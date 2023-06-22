#!/bin/bash

set -euo pipefail

input=$1
output=$2
threads=$3

prefix="${output%.*}"
raxml-ng \
    --all \
    --msa $input \
    --model GTR+G \
    --tree pars{1},rand{1} \
    --prefix $prefix \
    --bs-trees 5 \
    --threads $threads \
    --redo
cp $prefix.raxml.support $output
