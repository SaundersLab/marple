#!/bin/bash
# last update: 12/01/2024

run_cmd="snakemake --cores all --rerun-incomplete"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dev)
            run_cmd="/bin/bash"
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
		$run_cmd
        exit 0
	else
		echo "Environment marple-env not found."
        exit 1
	fi
else
    echo "Mamba not found."
    exit 1
fi
