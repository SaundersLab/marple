#!/bin/bash
# last update: 12/09/2024

run_cmd="snakemake --cores all --rerun-incomplete"
log_file="snakemake.log"

if command -v mamba &> /dev/null; then
    if mamba env list | grep -q marple-env; then
        mamba run -n marple-env snakemake --unlock &> /dev/null 
        mamba run -n marple-env $run_cmd > $log_file 2>&1 &
        pid=$!
        spin='-\|/'
        i=0
        while kill -0 $pid 2>/dev/null
        do
            i=$(( (i+1) %4 ))
            printf "\r${spin:$i:1}"
            sleep .1
        done
        wait $pid
        status=$?
        if [ $status -ne 0 ]; then
            echo "Command failed with exit status $status. Error log:"
            cat $log_file
            exit 1
        fi
        
        exit 0
    else
        echo "Conda environment marple-env not found."
        exit 1
    fi
else
    echo "Mamba not found."
    exit 1
fi
