#!/usr/bin/env bash

# Ensure correct environment
required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Parse arguments
dataset=""         # required: perfect or flip
scope="selected"   # default

while [[ $# -gt 0 ]]; do
  case $1 in
    --dataset)
      case ${2:-} in
        perfect|flip)
          dataset="simulated_pangenomes_${2}"
          ;;
        *)
          echo "Invalid dataset: ${2:-<missing>}"
          echo "Allowed values: perfect, flip"
          exit 1
          ;;
      esac
      shift 2
      ;;
    --scope)
      case ${2:-} in
        selected|all) scope="$2" ;;
        *)
          echo "Invalid scope: ${2:-<missing>}"
          echo "Allowed values: selected, all"
          exit 1
          ;;
      esac
      shift 2
      ;;
    *)
      echo "Usage: $0 --dataset <perfect|flip> [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

# Input directories
gpa_dir="real_pangenomes/gpa_matches${scope:+_${scope}}"
nodes_dir="real_pangenomes/coinfinder_runs_${scope}"

# Output base directory
out_base="${dataset}"
out_gpa_dir="${out_base}/gpa_matches${scope:+_${scope}}"

# Script path
script="scripts/pseudo-sim/sim_pan.py"

# Safety checks
for d in "$nodes_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found."
        exit 1
    fi
done

mkdir -p "$out_gpa_dir"

# Loop through each species
for nodes_file in "$nodes_dir"/*/coincident_nodes_all.tsv; do
    species_taxid=$(basename "$(dirname "$nodes_file")")
    gpa_file="${gpa_dir}/${species_taxid}_REDUCED.csv"

    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid} (scope=${scope})..."

    # Resolve paths
    nodes_file="$(realpath "$nodes_file")"
    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$out_base")"

    # Run simulation
    python3 "$script" \
        "$nodes_file" \
        "$gpa_file" \
        -o "${out_gpa_dir}/duplicates_${species_taxid}_REDUCED.csv" \
        --tab_out "${out_gpa_dir}/duplicates_${species_taxid}_REDUCED.tab" \
        --dup_mode "${dataset#simulated_pangenomes_}"

    echo "Finished ${species_taxid}"
    echo
done
