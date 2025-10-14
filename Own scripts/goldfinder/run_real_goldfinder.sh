#!/usr/bin/env bash

required_env="goldfinder"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Parse arguments
dataset=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --dataset)
      dataset="$2"
      shift 2
      ;;
    *)
      echo "Usage: $0 --dataset <real_pangenomes|simulated_pangenomes>"
      exit 1
      ;;
  esac
done

# Require dataset flag
if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real_pangenomes|simulated_pangenomes>"
    exit 1
fi

tree_dir="real_pangenomes/tree_matches"
gpa_dir="${dataset}/gpa_matches"
gold_dir="goldfinder/goldfinder"
out_base="${dataset}/goldfinder_runs"

# Safety: bail if directories aren’t found
for d in "$tree_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

mkdir -p "$out_base"

# Loop through each tree file in the tree_matches directory
for tree_file in "$tree_dir"/*.nwk; do
    
    filename=$(basename "$tree_file")
    basename_noext="${filename%.nwk}"
    species_taxid="${basename_noext#reduced_}"

    gpa_file=$(ls "${gpa_dir}"/*"${species_taxid}"_REDUCED.csv 2>/dev/null | head -n1)
    [ -f "$gpa_file" ] && echo "Found!" || echo "Missing!"

    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid}..."

    outdir="${out_base}/${species_taxid}"
    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "Skipping ${species_taxid}: $outdir already exists and is not empty."
        continue
    fi

    gpa_file="$(realpath "$gpa_file")"
    tree_file="$(realpath "$tree_file")"
    mkdir -p "$outdir"

    (
      NUMBA_NUM_THREADS=24 python3 "$gold_dir/goldfinder.py" \
        -i "$gpa_file" \
        -t "$tree_file" \
        -o "$outdir" \
        -g 50000 \
        -c both
    )

    echo "Finished ${species_taxid}"
done
