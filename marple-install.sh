#!/bin/bash
# last update: 16/05/2024

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        wd=/home/$USER/marple
        bashf=~/.bashrc
elif [[ "$OSTYPE" == "darwin"* ]]; then
        wd=/Users/$USER/marple
        bashf=~/.bash_profile
fi

cwd=$(pwd)

if [ $cwd != $wd ]; then
        echo "The marple directory should be in /home/$USER/marple. Move and try again"
        exit 1
fi

marple_func='\nfunction marple() {\npushd ~/marple > /dev/null 2>&1\nmamba activate marple-env \ncat .version.info \nbash marple.sh \nmamba deactivate \npopd > /dev/null 2>&1\n}'
transfer_pst_func='\nfunction transfer-pgt() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pgt/"$sample".fastq.gz\ndone\ndone\n}'
transfer_pgt_func='\nfunction transfer-pst() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pst/"$sample".fastq.gz\ndone\ndone\n}'
auspice_func='\nfunction marple-tree() {\nnpx kill-port 4000 > /dev/null 2>&1\nauspice develop --extend ~/marple/config/auspice/auspiceCustomisations/config.json --datasetDir ~/marple/results/auspice/ > /dev/null 2>&1 & xdg-open http://localhost:4000\n}'

# Create backup of .bashrc and add marple functions
if [[ ! -d ".marple-tmp" ]] && [[ $(grep -L transfer-pgt "$bashf") ]]; then
        mkdir .marple-tmp
        cp "$bashf" .marple-tmp/bashrc.bak
        echo -e "$marple_func" >> "$bashf"
        echo -e "$transfer_pst_func" >> "$bashf"
        echo -e "$transfer_pgt_func" >> "$bashf"
        echo -e "$auspice_func" >> "$bashf"
        echo -e '\nexport -f marple' >> "$bashf"
        echo -e 'export -f transfer-pst' >> "$bashf"
        echo -e 'export -f transfer-pgt' >> "$bashf"
        echo -e 'export -f marple-tree\n' >> "$bashf"
        source "$bashf"
elif [[ $(grep transfer-pgt "$bashf") ]]; then
        while true; do
        read -p "Do you want to update the functions for MARPLE? " yn
        case $yn in
        [Yy]* )
        mkdir -p .marple-tmp
        sed -i -e '/function marple() {/,/export -f marple-tree/d' "$bashf"
        mv "$bashf"-e .marple-tmp/bashrc.bak
        echo -e "$marple_func" >> "$bashf"
        echo -e "$transfer_pst_func" >> "$bashf"
        echo -e "$transfer_pgt_func" >> "$bashf"
        echo -e "$auspice_func" >> "$bashf"
        echo -e '\nexport -f marple' >> "$bashf"
        echo -e 'export -f transfer-pst' >> "$bashf"
        echo -e 'export -f transfer-pgt' >> "$bashf"
        echo -e 'export -f marple-tree\n' >> "$bashf"
        break;;
        [Nn]* ) break;;
        * ) echo "Please select an option [Yn] ";;
        esac
        done
fi

pckg="python=3.10 augur snakemake bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython"

update_software() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        sudo apt-get install bzip2
        sudo apt-get install gnumeric
        sudo apt install npm
        sudo npm install --global auspice
        sudo npm install react
        sudo npm install styled-components
        sudo npm install kill-port
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        brew install bzip2
        brew install gnumeric
        brew install npm
        npm install --global auspice
        npm install react
        npm install styled-components
        npm install kill-port
    fi
}

if command -v mamba &> /dev/null; then
        if mamba env list | grep -q marple-env; then
                while true; do
                read -p "Do you want to update the marple-env? " yn
                case $yn in
                [Yy]* )
                mamba env remove -n marple-env -y -q
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
                break;;
                [Nn]* ) break;;
                * ) echo "Please select an option [Yn] ";;
                esac
                done
        else
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
fi

while true; do
    read -p "Do you want to update the softwares used in MARPLE? " yn
    case $yn in
    [Yy]* ) update_software; break;;
    [Nn]* ) break;;
    * ) echo "Please select an option [Yn] ";;
    esac
done
