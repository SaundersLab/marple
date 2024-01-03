#!/bin/bash

set -euo pipefail

input=$1
output=$2

prefix="${output%.*}"
fasttree -gtr -nt < $1 > $2
