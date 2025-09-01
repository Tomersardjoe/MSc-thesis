#!/usr/bin/env bash

required_env="panforest"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

gpa_dir="real_pangenomes/gpa_matches"
pan_dir="panforest"
out_base="real_pangenomes/panforest_runs"

# Safety: bail if directory isn’t found
if [ ! -d "$gpa_dir" ]; then
    echo "Error: Directory '$gpa_dir' not found. Check the path."
    exit 1
fi

mkdir -p "$out_base"

# Loop through each GPA file in the gpa_matches directory
for gpa_file in "$gpa_dir"/*_REDUCED.csv; do
    
    # Extract species_taxid from filename
    filename=$(basename "$gpa_file")
    species_taxid="${filename%%_*}"

    echo "Starting ${species_taxid}..."

    # Output directory for this species
    outdir="${out_base}/${species_taxid}"
    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "Skipping ${species_taxid}: $outdir already exists and is not empty."
        continue
    fi
    mkdir -p "$outdir"

    # Use absolute file paths
    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    # Step 1: Process matrix
    python3 $pan_dir/process_matrix.py \
        -i "$gpa_file" \
        -o "${outdir}/collapsed_matrix.csv" \
        -d "$outdir"

    # Step 2: Run PanForest
    python3 $pan_dir/PanForest.py \
        -n 500 \
        -d 8 \
        -m "${outdir}/collapsed_matrix.csv" \
        -pres 1 \
        -abs 1 \
        -o "$outdir"

    echo "Finished ${species_taxid}"
done
