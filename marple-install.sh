#!/bin/bash
# last update: 02/04/2024

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

marple_func='\nfunction marple() {\npushd ~/marple \nmamba activate marple-env \ncat .version.info \nbash marple.sh $1 \nmamba deactivate \npopd\n}'
transfer_pst_func='\nfunction transfer-pgt() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pgt/"$sample".fastq.gz\ndone\ndone\n}'
transfer_pgt_func='\nfunction transfer-pst() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > ~/marple/reads/pst/"$sample".fastq.gz\ndone\ndone\n}'
auspice_func='\nfunction view-marple-tree() {\nauspice view --datasetDir ~/marple/results/auspice/ & xdg-open http://localhost:4000\n}'

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
        echo -e 'export -f view-marple-tree\n' >> "$bashf"
        source "$bashf"
elif [[ $(grep transfer-pgt "$bashf") ]]; then
        while true; do
        read -p "Do you want to update the functions for MARPLE? " yn
        case $yn in
        [Yy]* )
        mkdir -p .marple-tmp
        sed -i-e '/function marple() {/,/export -f view-marple-tree/d' "$bashf";
        mv "$bashf"-e .marple-tmp/bashrc.bak
        echo -e "$marple_func" >> "$bashf";
        echo -e "$transfer_pst_func" >> "$bashf";
        echo -e "$transfer_pgt_func" >> "$bashf";
        echo -e "$auspice_func" >> "$bashf"
        echo -e '\nexport -f marple' >> "$bashf"
        echo -e 'export -f transfer-pst' >> "$bashf"
        echo -e 'export -f transfer-pgt' >> "$bashf"
        echo -e 'export -f view-marple-tree\n' >> "$bashf"
        break;;
        [Nn]* ) break;;
        * ) echo "Please select an option [Yn] ";;
        esac
        done
fi

pckg="python=3.10 augur snakemake bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython"

if command -v mamba &> /dev/null; then
        if mamba env list | grep -q marple-env; then
                while true; do
                read -p "Do you want to update the marple-env? " yn
                case $yn in
                [Yy]* )
                mamba env remove -n marple-env -y -q
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
                break;;
                [Nn]* ) while true; do
                read -p "Do you want to update the softwares used in MARPLE? " yn
                case $yn in
                [Yy]* ) 
                if [[ "$OSTYPE" == "linux-gnu"* ]]; then
                sudo apt-get install bzip2
                sudo apt-get install gnumeric
                sudo apt install npm
                sudo npm install --global auspice
                elif [[ "$OSTYPE" == "darwin"* ]]; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
                brew install bzip2
                brew install gnumeric
                brew install npm
                npm install --global auspice
                fi 
                break;;
                [Nn]* ) break;;
                * ) echo "Please select an option [Yn] ";;
                esac
                done
                break;;
                * ) echo "Please select an option [Yn] ";;
                esac
                done
        else
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
        fi
elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "darwin"* ]]; then
        if [[ "$OSTYPE" == "linux-gnu"* ]] && [[ ! -d ~/micromamba ]]; then
                sudo apt-get install bzip2
                sudo apt-get install gnumeric
                sudo apt install npm
                sudo npm install --global auspice
                wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
                ./bin/micromamba shell init -s bash -p ~/micromamba
                source "$bashf"
                sleep 1
                micromamba create -n base -c conda-forge -y
                micromamba activate base
        elif [[ "$OSTYPE" == "darwin"* ]] && ! command -v micromamba &> /dev/null; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
                brew install bzip2
                brew install micromamba
                brew install gnumeric
                brew install npm
                npm install --global auspice
                curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/mac | bash
                eval "$(micromamba shell hook --shell bash)"
                micromamba activate
        fi
        micromamba install mamba -c conda-forge -y
        mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
        mamba init
else
        echo "Unsupported operating system"
        exit 1
fi
