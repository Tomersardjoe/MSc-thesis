#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

nodes_dir="real_pangenomes/coinfinder_runs"
gpa_dir="real_pangenomes/gpa_matches"
script="scripts/pseudo-sim/sim_pan.py"
out_base="simulated_pangenomes"

# Safety checks
for d in "$nodes_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found."
        exit 1
    fi
done

mkdir -p "$out_base"

# Loop through each nodes_all.tsv file
for nodes_file in "$nodes_dir"/*/coincident_nodes_all.tsv; do
    # Extract species_taxid from parent directory name
    species_taxid=$(basename "$(dirname "$nodes_file")")

    # Find matching GPA file
    gpa_file="$gpa_dir/${species_taxid}_REDUCED.csv"
    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid}..."

    # Output directory
    outdir="${out_base}"
    mkdir -p "$outdir/gpa_matches"

    # Absolute paths
    nodes_file="$(realpath "$nodes_file")"
    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    # Run simulation script
    python3 "$script" \
        "$nodes_file" \
        "$gpa_file" \
        -o "${outdir}/gpa_matches/duplicates_${species_taxid}_REDUCED.csv" \
        --tab_out "${outdir}/gpa_matches/duplicates_${species_taxid}_REDUCED.tab"

    echo "Finished ${species_taxid}"
done
