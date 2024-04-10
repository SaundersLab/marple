#!/bin/bash
# last update: 10/04/2024

run_cmd="nohup snakemake --quiet --cores all --rerun-incomplete"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev)
            run_cmd="/bin/bash; conda activate marple-env"
            shift
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

if command -v mamba &> /dev/null; then
    if mamba env list | grep -q marple-env; then
    	conda activate marple-env > /dev/null 2>&1
		$run_cmd 2> /dev/null &
        pid=$!
        spin='-\|/'
        i=0
        while kill -0 $pid 2>/dev/null
        do
            i=$(( (i+1) %4 ))
            printf "\r${spin:$i:1}"
            sleep .1
        done 
        mv nohup.out logs/snakemake.out 
        exit 0
	else
		echo "Environment marple-env not found."
        exit 1
	fi
else
    echo "Mamba not found."
    exit 1
fi
