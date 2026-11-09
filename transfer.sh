#!/bin/bash
# Concatenates the fastq files of a MinKNOW experiment into reads/<organism>/<sample>.fastq.gz.
# Run through the transfer-pst and transfer-pgt functions that marple-install.sh adds to .bashrc.

marple_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
minknow_data="${MINKNOW_DATA:-/var/lib/minknow/data}"

die() {
    echo "Error: $*" >&2
    exit 1
}

usage() {
    cat >&2 <<EOF
Usage:
  transfer-$organism <name of experiment> --samplesheet <samplesheet.csv>
  transfer-$organism <name of experiment> <barcodeXX>=<name of sample> <barcodeYY>=<name of sample> ...

The samplesheet needs a header row with the columns barcode and sample:
  barcode,sample
  barcode01,field1
  barcode02,field2
EOF
    exit 1
}

# Prints one barcode=sample line per row. Handles the byte order mark and
# Windows line endings that Excel adds when it saves a csv.
samplesheet_to_pairs() {
    LC_ALL=C awk -F',' '
        function clean(field) {
            gsub(/^[ \t\r"]+|[ \t\r"]+$/, "", field)
            return field
        }
        NR == 1 {
            if (substr($0, 1, 3) == "\357\273\277") $0 = substr($0, 4)
            for (i = 1; i <= NF; i++) {
                column = tolower(clean($i))
                if (column == "barcode") barcode_col = i
                if (column == "sample") sample_col = i
            }
            if (!barcode_col || !sample_col) exit
            next
        }
        {
            barcode = clean($barcode_col)
            sample = clean($sample_col)
            if (barcode != "" || sample != "") print barcode "=" sample
        }
        END {
            if (!barcode_col || !sample_col) {
                print "Error: the samplesheet needs a header row with the columns barcode and sample" > "/dev/stderr"
                exit 1
            }
        }
    ' "$1"
}

# Prints every directory of reads for the experiment, one per line. Live
# basecalling writes to fastq_pass, and basecalling after a run writes to
# basecalling/pass.
find_read_dirs() {
    find -H "$minknow_data" -type d -name "$experiment" \
        -exec find {} -type d \( -name fastq_pass -o -path '*/basecalling/pass' \) \; 2>/dev/null | sort -u
}

transfer_sample() {
    local barcode=$1 sample=$2
    local dest="$marple_dir/reads/$organism/$sample.fastq.gz"
    set -- "$read_dir/$barcode"/*.fastq.gz
    if [[ ! -e "$1" ]]; then
        echo "  FAILED $barcode -> $sample: no .fastq.gz files in $read_dir/$barcode"
        return 1
    fi
    if cat "$@" > "$dest.tmp" && mv -f "$dest.tmp" "$dest"; then
        echo "  $barcode -> $sample ($# files)"
    else
        rm -f "$dest.tmp"
        echo "  FAILED $barcode -> $sample: could not write $dest"
        return 1
    fi
}

organism=$1
experiment=$2
if [[ "$experiment" == "" || "$experiment" == -* || $# -lt 3 ]]; then
    usage
fi
shift 2

pairs=()
if [[ "$1" == "--samplesheet" || "$1" == "-s" ]]; then
    [[ $# -eq 2 ]] || usage
    [[ -r "$2" ]] || die "cannot read samplesheet: $2"
    samplesheet_pairs=$(samplesheet_to_pairs "$2") || exit 1
    [[ -n "$samplesheet_pairs" ]] || die "no samples found in $2"
    while IFS= read -r pair; do
        pairs+=("$pair")
    done <<< "$samplesheet_pairs"
else
    pairs=("$@")
fi

# Check every pair before transferring anything
invalid=0
barcodes=()
samples=()
for pair in "${pairs[@]}"; do
    barcode=${pair%%=*}
    sample=${pair#*=}
    if [[ "$pair" != *=* || ! "$barcode" =~ ^[A-Za-z0-9._-]+$ || ! "$sample" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "Error: invalid barcode and sample '$pair'. Names can only contain letters, numbers, '.', '_' and '-'." >&2
        invalid=1
    fi
    barcodes+=("$barcode")
    samples+=("$sample")
done
duplicates=$(printf '%s\n' "${barcodes[@]}" | sort | uniq -d)
if [[ -n "$duplicates" ]]; then
    echo "Error: barcodes listed more than once: $(echo $duplicates)" >&2
    invalid=1
fi
duplicates=$(printf '%s\n' "${samples[@]}" | sort | uniq -d)
if [[ -n "$duplicates" ]]; then
    echo "Error: sample names listed more than once: $(echo $duplicates)" >&2
    invalid=1
fi
[[ "$invalid" -eq 0 ]] || exit 1

# Use the most recent directory of reads
read_dir=""
read_dir_count=0
while IFS= read -r dir; do
    read_dir_count=$((read_dir_count + 1))
    if [[ -z "$read_dir" || "$dir" -nt "$read_dir" ]]; then
        read_dir=$dir
    fi
done < <(find_read_dirs)
if [[ -z "$read_dir" ]]; then
    die "no fastq_pass or basecalling/pass directory found for experiment '$experiment' in $minknow_data"
fi
if [[ "$read_dir_count" -gt 1 ]]; then
    echo "Found $read_dir_count directories of reads for experiment '$experiment'. Only the most recent is used:"
    find_read_dirs | sed 's/^/  /'
fi
echo "Transferring reads from $read_dir"

mkdir -p "$marple_dir/reads/$organism" || die "cannot create $marple_dir/reads/$organism"
transferred=0
for i in "${!barcodes[@]}"; do
    transfer_sample "${barcodes[$i]}" "${samples[$i]}" && transferred=$((transferred + 1))
done

echo "Transferred $transferred of ${#barcodes[@]} samples to $marple_dir/reads/$organism"
[[ "$transferred" -eq "${#barcodes[@]}" ]]
