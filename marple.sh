#!/bin/bash
# last update: 09/09/2024

run_cmd="nohup snakemake --cores all --rerun-incomplete"

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
        exit 0
	else
		echo "Environment marple-env not found."
        exit 1
	fi
else
    echo "Mamba not found."
    exit 1
fi
