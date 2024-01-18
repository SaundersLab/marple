#!/bin/bash
# last update: 12/01/2024

wd=/home/$USER/marple
cwd=$(pwd)

if [ $cwd != $wd ]; then
        echo "The marple directory should be in /home/$USER/marple. Move and try again"
        :
fi

pckg="snakemake bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython"

if ! grep transfer-pgt ~/.bashrc &> /dev/null; then
        echo -e '\nfunction transfer-pgt () {\ncat /var/lib/minknow/data/reads/$1/*/*/basecalling/pass/$2/*.fastq.gz > /home/$USER/marple/reads/pgt/$3.fastq.gz\n}' >> ~/.bashrc
        echo -e '\nfunction transfer-pst () {\ncat /var/lib/minknow/data/reads/$1/*/*/basecalling/pass/$2/*.fastq.gz > /home/$USER/marple/reads/pst/$3.fastq.gz\n}' >> ~/.bashrc
        echo -e '\nfunction marple() {\npushd /home/$USER/marple \nmamba activate marple-env \nbash marple.sh $1 \nmamba deactivate \npopd\n}' >> ~/.bashrc
        source ~/.bashrc
fi

if command -v mamba &> /dev/null; then
if mamba env list | grep -q marple-env; then
                :
        else
                mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
        fi
elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "darwin"* ]]; then
        echo "$OSTYPE"
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
                sudo apt-get install bzip2
        elif [[ "$OSTYPE" == "darwin"* ]]; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
                brew install bzip2
        fi
        wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
        ./bin/micromamba shell init -s bash -p ~/micromamba
        source ~/.bashrc
        sleep 1
        micromamba create -n base -c conda-forge -y
        micromamba activate base
        micromamba install mamba -c conda-forge -y
        mamba create -n marple-env -y -c bioconda -c conda-forge $pckg
        mamba init
else
        echo "Unsupported operating system"
        exit 1
fi
