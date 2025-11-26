#!/bin/bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

base_outdir="${PARENT_DIR}/real_pangenomes"
mkdir -p "$base_outdir"

gpa_dir="${base_outdir}/gpa_matches_all"
output_csv="${base_outdir}/alpha_fluidity_all.csv"

# Write header if file doesn’t exist
if [ ! -f "$output_csv" ]; then
    echo "species_taxid,filename,fluidity_calc,openness" > "$output_csv"
fi

# Safety: bail if directories aren’t found
if [ ! -d "$gpa_dir" ]; then
    echo "Error: Directory '$gpa_dir' not found. Check the path."
    exit 1
fi

for file in "$gpa_dir"/*_REDUCED.csv; do
    filename=$(basename "$file")
    species_taxid=$(echo "$filename" | cut -d'_' -f1)

    # Skip if already processed
    pattern=$(printf '%s' "$species_taxid" | sed -e 's/[]\/$*.^|[]/\\&/g')
    if sed 's/\r$//' "$output_csv" | grep -Eq "^$pattern,"; then
        echo "Skipping $filename..."
        continue
    fi

    echo "Processing $filename..."
    result=$(Rscript "$SCRIPT_DIR/fluidity_calc.R" "$file")

    fluidity=$(echo "$result" | grep "Genomic fluidity" | sed -E 's/.*: *([0-9.]+) *$/\1/')
    openness=$(echo "$result" | grep "Pangenome openness" | sed -E 's/.*: *([0-9.]+) *$/\1/')

    echo "$species_taxid,$filename,$fluidity,$openness" >> "$output_csv"
done

echo "All results written to $output_csv"
