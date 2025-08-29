#!/usr/bin/env bash

required_env="goldfinder"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

tree_dir="real_pangenomes/tree_matches"
gpa_dir="real_pangenomes/gpa_matches"
gold_dir="goldfinder/goldfinder
out_base="real_pangenomes/goldfinder_runs"

# Safety: bail if directories aren’t found
for d in "$tree_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

mkdir -p "$out_base"

# Loop through each tree file in the tree_matches directory
for tree_file in "$tree_dir"/*_red_tree_converted.nwk; do
    
    # Extract species_taxid from filename
    filename=$(basename "$tree_file")
    species_taxid="${filename%%_*}"

    # Build path to corresponding GPA file
    gpa_file="${gpa_dir}/${species_taxid}_REDUCED.csv"
    [ -f "$gpa_file" ] && echo "Found!" || echo "Missing!"

    # Skip if GPA file doesn't exist
    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid}..."

    # Output directory for this species
    outdir="${out_base}/${species_taxid}"
    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "Error: $outdir already exists and is not empty. Stopping."
        exit 1
    fi

    # Use absolute file paths
    gpa_file="$(realpath "$gpa_file")"
    tree_file="$(realpath "$tree_file")"
    
    mkdir -p "$outdir"

    # Run Goldfinder
    (
      NUMBA_NUM_THREADS=24 python3 $gold_dir/goldfinder.py \
        -i "$gpa_file" \
        -t "$tree_file" \
        -o "$outdir" \
        -g 50000 \
        -c both
    )

    echo "Finished ${species_taxid}"
done
