#!/bin/bash

# Absolute paths
base_dir="$(pwd)"
script_path="$(realpath ../simulation/fluidity_calc.R)"
output_csv="$base_dir/alpha_fluidity_all.csv"

if [ ! -f "$output_csv" ]; then
    echo "filename,fluidity,openness" > "$output_csv"
fi

for file in "$base_dir"/gpa_matches/*_REDUCED.csv; do
    filename=$(basename "$file")

    # Skip if this filename and both numeric values are present already
    pattern=$(printf '%s' "$filename" | sed -e 's/[]\/$*.^|[]/\\&/g')
    if sed 's/\r$//' "$output_csv" | grep -Eq "^$pattern,[0-9.]+,[0-9.]+$"; then
    echo "Skipping..."
    continue
    fi

    echo "Processing $filename..."
    result=$(Rscript "$script_path" "$file")

    # Prevent grepping progress lines
    fluidity=$(echo "$result" | grep "Genomic fluidity" | sed -E 's/.*: *([0-9.]+) *$/\1/')
    openness=$(echo "$result" | grep "Pangenome openness" | sed -E 's/.*: *([0-9.]+) *$/\1/')

    echo "$filename,$fluidity,$openness" >> "$output_csv"
done

echo "All results written to $output_csv"
