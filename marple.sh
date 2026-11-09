#!/bin/bash
# last update: 09/11/2026

run_cmd="snakemake --cores all --rerun-incomplete"
log_file="snakemake.log"
status_file="results/pipeline_status.txt"

# Keep running if the terminal closes. Snakemake inherits this setting.
trap '' HUP

# Records the outcome in a file and a desktop notification, which are still
# there if the terminal has closed.
report_outcome() {
    echo "$1 ($(date))" > "$status_file"
    if command -v notify-send &> /dev/null; then
        notify-send -u critical "MARPLE" "$1" 2> /dev/null
    fi
}

if command -v mamba &> /dev/null; then
    if mamba env list | grep -q marple-env; then
        mkdir -p results
        echo "Pipeline started and not yet finished ($(date))" > "$status_file"
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
        printf "\r \r"
        if [ $status -ne 0 ]; then
            report_outcome "Pipeline failed with exit status $status. Error log: $PWD/$log_file"
            echo "Command failed with exit status $status. Error log:"
            cat $log_file
            exit 1
        fi

        report_outcome "Pipeline completed successfully. Results can be found in $PWD/results"
        echo "Pipeline completed successfully. Results can be found in $PWD/results"
        for report in results/*/report/*.multiqc.html; do
            [ -e "$report" ] && echo "  Report: $PWD/$report"
        done
        exit 0
    else
        echo "Conda environment marple-env not found."
        exit 1
    fi
else
    echo "Mamba not found."
    exit 1
fi
