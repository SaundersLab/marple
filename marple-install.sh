#!/bin/bash
# last update: 09/11/2026

# This script is normally sourced, and exit would close the user's terminal.
# Every early stop therefore uses: return 1 2>/dev/null || exit 1

wd="$HOME/marple"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    bashf="$HOME/.bashrc"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    bashf="$HOME/.bash_profile"
else
    echo "Unsupported OS type: $OSTYPE"
    return 1 2>/dev/null || exit 1
fi

# Check marple is in the correct directory
cwd=$(pwd)
if [[ "$cwd" != "$wd" ]]; then
    echo "The marple directory should be in $wd, but you are in $cwd. Move it there and try again."
    return 1 2>/dev/null || exit 1
fi

# Check dependencies before changing anything
missing_deps=0
if ! command -v curl &> /dev/null; then
    missing_deps=1
    echo "curl not found. Install it, then run this script again:"
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "    sudo apt-get install curl"
    else
        echo "    brew install curl"
    fi
    echo "    Other options: https://curl.se/download.html"
fi
if ! command -v mamba &> /dev/null; then
    missing_deps=1
    echo "mamba not found. Install Miniforge, which includes mamba, then open a new terminal and run this script again:"
    echo "    https://github.com/conda-forge/miniforge#install"
fi
if [[ "$missing_deps" -ne 0 ]]; then
    echo "MARPLE was not installed."
    return 1 2>/dev/null || exit 1
fi

marple_func='\nfunction marple() {\npushd ~/marple > /dev/null 2>&1 \neval "$(mamba shell hook --shell bash)" \ncat .version.info \nbash marple.sh \npopd > /dev/null 2>&1\n}'
transfer_pst_func='
function transfer-pst() {
  experiment=$1
  shift
  for arg in "$@"; do
    IFS="=" read -r barcode sample <<< "$arg"
    rm -f ~/marple/reads/pst/"$sample".fastq.gz
    # Find the most recent basecalling directory
    newest_dir=$(find /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | \
                 xargs -I {} find "{}" -type d -name basecalling | \
                 xargs -d "\n" ls -dt 2>/dev/null | head -n 1)
    if [ -n "$newest_dir" ]; then
      cat "$newest_dir/pass/$barcode"/*.fastq.gz >> ~/marple/reads/pst/"$sample".fastq.gz 2>/dev/null
    else
      echo "No basecalling directory found for experiment: $experiment"
    fi
  done
}'
transfer_pgt_func='
function transfer-pgt() {
  experiment=$1
  shift
  for arg in "$@"; do
    IFS="=" read -r barcode sample <<< "$arg"
    rm -f ~/marple/reads/pgt/"$sample".fastq.gz
    # Find the most recent basecalling directory
    newest_dir=$(find /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | \
                 xargs -I {} find "{}" -type d -name basecalling | \
                 xargs -d "\n" ls -dt 2>/dev/null | head -n 1)
    if [ -n "$newest_dir" ]; then
      cat "$newest_dir/pass/$barcode"/*.fastq.gz >> ~/marple/reads/pgt/"$sample".fastq.gz 2>/dev/null
    else
      echo "No basecalling directory found for experiment: $experiment"
    fi
  done
}'

# Create backup of .bashrc and add marple functions
if [[ ! -d ".marple-tmp" ]] && ! grep -q 'transfer-pgt' "$bashf"; then
    mkdir .marple-tmp
    cp "$bashf" .marple-tmp/bashrc.bak
    {
        echo -e "$marple_func"
        echo -e "$transfer_pst_func"
        echo -e "$transfer_pgt_func"
        echo -e '\nexport -f marple'
        echo -e 'export -f transfer-pst'
        echo -e 'export -f transfer-pgt'
    } >> "$bashf"
    source "$bashf"
elif grep -q 'transfer-pgt' "$bashf"; then
    while true; do
        read -p "Do you want to update the functions for MARPLE? (y/n) " yn
        case $yn in
            [Yy]* )
                mkdir -p .marple-tmp
                sed -i -e '/function marple() {/,/export -f transfer-pgt/d' "$bashf"
                [ -f "$bashf"-e ] && mv "$bashf"-e .marple-tmp/bashrc.bak
                {
                    echo -e "$marple_func"
                    echo -e "$transfer_pst_func"
                    echo -e "$transfer_pgt_func"
                    echo -e '\nexport -f marple'
                    echo -e 'export -f transfer-pst'
                    echo -e 'export -f transfer-pgt'
                } >> "$bashf"
                source "$bashf"
                break;;
            [Nn]* ) break;;
            * ) echo "Please select an option [y/n] ";;
        esac
    done
fi

update_software() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        sudo apt-get update
        sudo apt-get install -y bzip2 gnumeric
        # clustalw2 installation
        sudo apt-get install -y build-essential
        curl -fsSL -O https://ftp.ebi.ac.uk/pub/software/clustalw2/2.1/clustalw-2.1.tar.gz
        (tar xzf clustalw-2.1.tar.gz && cd clustalw-2.1/ && ./configure && make && sudo make install)
        rm -rf clustalw-2.1 clustalw-2.1.tar.gz
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew install bzip2 gnumeric
    fi
}

if mamba env list | grep -q marple-env; then
    while true; do
        read -p "Do you want to update the marple-env? (y/n) " yn
        case $yn in
            [Yy]* )
                mamba env remove -n marple-env -y -q
                mamba env create -f config/env.yml
                break;;
            [Nn]* ) break;;
            * ) echo "Please select an option [y/n] ";;
        esac
    done
else
    mamba env create -f config/env.yml
fi

while true; do
    read -p "Do you want to install/update the softwares used in MARPLE? (y/n) " yn
    case $yn in
        [Yy]* ) update_software; break;;
        [Nn]* ) break;;
        * ) echo "Please select an option [y/n] ";;
    esac
done

# Report on what is installed, whichever options were chosen above
install_problems=()
grep -q 'transfer-pgt' "$bashf" || install_problems+=("MARPLE functions are missing from $bashf")
mamba run -n marple-env snakemake --version &> /dev/null || install_problems+=("snakemake does not run in the marple-env environment")
for tool in clustalw2 bzip2 ssconvert; do
    command -v "$tool" &> /dev/null || install_problems+=("$tool not found")
done

echo
if [[ ${#install_problems[@]} -eq 0 ]]; then
    echo "MARPLE installed successfully. Open a new terminal, or run: source $bashf"
else
    echo "MARPLE installation is incomplete:"
    printf '  - %s\n' "${install_problems[@]}"
    echo "Check the output above for errors, then run this script again."
    return 1 2>/dev/null || exit 1
fi
