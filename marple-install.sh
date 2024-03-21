#!/bin/bash
# last update: 21/03/2024

wd=/home/$USER/marple
cwd=$(pwd)

if [ $cwd != $wd ]; then
        echo "The marple directory should be in /home/$USER/marple. Move and try again"
        exit 1
fi

# These old transfer-p[g/s]t may be useful for future GUI implementation
# echo -e '\nfunction transfer-pgt () {\ncat /var/lib/minknow/data/reads/$1/*/*/basecalling/pass/$2/*.fastq.gz > /home/$USER/marple/reads/pgt/$3.fastq.gz\n}' >> ~/.bashrc
# echo -e '\nfunction transfer-pst () {\ncat /var/lib/minknow/data/reads/$1/*/*/basecalling/pass/$2/*.fastq.gz > /home/$USER/marple/reads/pst/$3.fastq.gz\n}' >> ~/.bashrc

marple_func='\nfunction marple() {\npushd /home/$USER/marple \nmamba activate marple-env \ncat .version.info \nbash marple.sh $1 \nmamba deactivate \npopd\n}'
transfer_pst_func='\nfunction transfer-pgt() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > /home/$USER/marple/reads/pgt/"$sample".fastq.gz\ndone\ndone\n}'
transfer_pgt_func='\nfunction transfer-pst() {\nexperiment=$1\nshift\nfor arg in "$@";do\nIFS="=" read -r barcode sample<<<"$arg"\nfind /var/lib/minknow/data/* -type d -name "$experiment" 2>/dev/null | xargs -I {} find '{}' -type d -name basecalling | while read dir; do\ncat $dir/pass/"$barcode"/*.fastq.gz > /home/$USER/marple/reads/pst/"$sample".fastq.gz\ndone\ndone\n}'

# Create backup of .bashrc and add marple functions
if [[ ! -d ".marple-tmp" ]] && [[ $(grep -L transfer-pgt ~/.bashrc) ]]; then
        mkdir .marple-tmp
        cp ~/.bashrc .marple-tmp/bashrc.bak
        echo -e "$marple_func" >> ~/.bashrc
        echo -e "$transfer_pst_func" >> ~/.bashrc
        echo -e "$transfer_pgt_func" >> ~/.bashrc
        echo -e '\nexport -f marple' >> ~/.bashrc
        echo -e 'export -f transfer-pst' >> ~/.bashrc
        echo -e 'export -f transfer-pgt\n' >> ~/.bashrc
        source ~/.bashrc
elif [[ $(grep transfer-pgt ~/.bashrc) ]]; then
        echo "Do you want to update the functions for MARPLE?"
        select yn in "Yes" "No"; do
        case $yn in
        Yes )
        mkdir -p .marple-tmp
        sed -i-e '/function marple() {/,/export -f transfer-pgt/d' ~/.bashrc;
        mv ~/.bashrc-e .marple-tmp/bashrc.bak
        echo -e "$marple_func" >> ~/.bashrc;
        echo -e "$transfer_pst_func" >> ~/.bashrc;
        echo -e "$transfer_pgt_func" >> ~/.bashrc;
        echo -e '\nexport -f marple' >> ~/.bashrc
        echo -e 'export -f transfer-pst' >> ~/.bashrc
        echo -e 'export -f transfer-pgt\n' >> ~/.bashrc
        break;;
        No ) break;;
        esac
        done
fi

pckg="snakemake bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython"

if command -v mamba &> /dev/null; then
        if mamba env list | grep -q marple-env; then
                :
        else
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
        fi
elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "darwin"* ]]; then
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
                sudo apt-get install bzip2
                wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
                ./bin/micromamba shell init -s bash -p ~/micromamba
                source ~/.bashrc
                sleep 1
                micromamba create -n base -c conda-forge -y
                micromamba activate base
        elif [[ "$OSTYPE" == "darwin"* ]]; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
                brew install bzip2
                brew install micromamba
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
