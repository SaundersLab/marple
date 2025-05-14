#!/bin/bash
# last update: 14/05/2025

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    wd="/home/$USER/marple"
    bashf="$HOME/.bashrc"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wd="/Users/$USER/marple"
    bashf="$HOME/.bash_profile"
else
    echo "Unsupported OS type: $OSTYPE"
    exit 1
fi

# Check marple is in the correct directory
cwd=$(pwd)
if [[ "$cwd" != "$wd" ]]; then
    echo "The marple directory should be in $wd. Move and try again"
    exit 1
fi

marple_func='\nfunction marple() {\npushd ~/marple > /dev/null 2>&1 \neval "$(mamba shell hook --shell bash)" \nmamba activate marple-env \ncat .version.info \nbash marple.sh \nmamba deactivate \npopd > /dev/null 2>&1\n}'
transfer_pst_func='\nfunction transfer-pst() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pst/"$sample".fastq.gz\ndone\ndone\n}'
transfer_pgt_func='\nfunction transfer-pgt() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pgt/"$sample".fastq.gz\ndone\ndone\n}'

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
        sudo apt-get install -y gpp g++ c++ kcc fcc gpp
        wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz
        tar xvzf clustalw-2.1.tar.gz
        cd clustalw-2.1/
        ./configure
        make
        sudo make install
        cd ../; rm -rf clustalw*
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew install bzip2 gnumeric
    fi
}

if command -v mamba &> /dev/null; then
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
else
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate
        micromamba install mamba -c conda-forge -y
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew install micromamba
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate
        micromamba install mamba -c conda-forge -y
    fi
    mamba shell init --shell bash
    source "$bashf"
    eval "$(mamba shell hook --shell bash)"
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
