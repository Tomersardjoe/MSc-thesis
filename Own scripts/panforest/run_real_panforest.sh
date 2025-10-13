#!/usr/bin/env bash

required_env="panforest"
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

gpa_dir="${dataset}/gpa_matches"
pan_dir="panforest"
out_base="${dataset}/panforest_runs"

# Safety: bail if directory isn’t found
if [ ! -d "$gpa_dir" ]; then
    echo "Error: Directory '$gpa_dir' not found. Check the path."
    exit 1
fi

mkdir -p "$out_base"

# Loop through each GPA file in the gpa_matches directory
for gpa_file in "$gpa_dir"/*_REDUCED.csv; do
    
    filename=$(basename "$gpa_file")
    # Remove the suffix
    base="${filename%_REDUCED.csv}"
    # Take the part after the last underscore
    species_taxid="${base##*_}"

    echo "Starting ${species_taxid}..."

    outdir="${out_base}/${species_taxid}"
    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "Skipping ${species_taxid}: $outdir already exists and is not empty."
        continue
    fi
    mkdir -p "$outdir"

    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    # Step 1: Process matrix
    python3 "$pan_dir/process_matrix.py" \
        -i "$gpa_file" \
        -o "${outdir}/collapsed_matrix.csv" \
        -d "$outdir"

    # Step 2: Run PanForest
    python3 "$pan_dir/PanForest.py" \
        -n 1000 \
        -d 16 \
        -m "${outdir}/collapsed_matrix.csv" \
        -pres 1 \
        -abs 1 \
        -o "$outdir"

    echo "Finished ${species_taxid}"
done
