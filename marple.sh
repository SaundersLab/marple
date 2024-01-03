#!/bin/bash
# last update: 03/01/2024

cwd=$(pwd)

run_cmd="docker run -v $cwd:/workflow -w /workflow marple snakemake --cores all"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev)
            run_cmd="docker run -v $cwd:/workflow -w /workflow -it marple /bin/bash"
            shift
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

if command -v docker &> /dev/null; then
    echo "Skipping Docker installation"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "Detected Linux"
    sudo apt update
    sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
    echo "deb [signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt update
    sudo apt install -y docker-ce docker-ce-cli containerd.io
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install --cask docker
else
    echo "Unsupported operating system"
    exit 1
fi

if command -v git &> /dev/null; then
    echo "Skipping git installation"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "Detected Linux"
    sudo apt update
    sudo apt install git 
elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS"
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install git
else
    echo "Unsupported operating system"
    exit 1
fi

docker --version
git --version

# Create a docker img of snakemake if it doesn't exist

if docker images | grep marple > /dev/null; then
	$run_cmd
else
	git clone https://github.com/snakemake/snakemake.git
	cd snakemake/
	if [[  "$OSTYPE" == "linux-gnu"* ]]; then
		sed -i 's/apptainer/apptainer -c bioconda bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython/g' Dockerfile
		sudo groupadd -f docker
		sudo usermod -aG docker $USER
    		sg docker -c '
       		docker build -t marple ./
        	cd ..
        	rm -rf snakemake
        	$run_cmd'
	elif [[ "$OSTYPE" == "darwin"* ]]; then
		sed -i -tmp 's/apptainer/apptainer -c bioconda bwa star samtools nanoq fastqc gffread multiqc fasttree openpyxl matplotlib biopython/g' Dockerfile
		rm Dockerfile-tmp
		docker build -t marple ./
		cd ..
		rm -rf snakemake
		$run_cmd
	else
		echo "Unsupported operating system"
	fi
fi

