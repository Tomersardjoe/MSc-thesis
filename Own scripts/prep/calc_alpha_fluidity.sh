#!/bin/bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Locate directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# Top-level output directory under repo root
base_outdir="${PARENT_DIR}/real_pangenomes"
mkdir -p "$base_outdir"

gpa_dir="$(realpath "$base_outdir/gpa_matches")"
output_csv="$base_outdir/alpha_fluidity_all.csv"

if [ ! -f "$output_csv" ]; then
    echo "filename,fluidity,openness" > "$output_csv"
fi

# Safety: bail if directories aren’t found
for d in "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

for file in "$base_outdir"/gpa_matches_all/*_REDUCED.csv; do
    filename=$(basename "$file")

    # Skip if this filename and both numeric values are present already
    pattern=$(printf '%s' "$filename" | sed -e 's/[]\/$*.^|[]/\\&/g')
    if sed 's/\r$//' "$output_csv" | grep -Eq "^$pattern,[0-9.]+,[0-9.]+$"; then
    echo "Skipping..."
    continue
    fi

    echo "Processing $filename..."
    result=$(Rscript "$SCRIPT_DIR/fluidity_calc.R" "$file")

    # Prevent grepping progress lines
    fluidity=$(echo "$result" | grep "Genomic fluidity" | sed -E 's/.*: *([0-9.]+) *$/\1/')
    openness=$(echo "$result" | grep "Pangenome openness" | sed -E 's/.*: *([0-9.]+) *$/\1/')

    echo "$filename,$fluidity,$openness" >> "$output_csv"
done

echo "All results written to $output_csv"
