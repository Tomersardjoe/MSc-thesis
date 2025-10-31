#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments
dataset=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --dataset)
      case ${2:-} in
        real)    dataset="real_pangenomes" ;;
        perfect) dataset="simulated_pangenomes_perfect" ;;
        flip)    dataset="simulated_pangenomes_flip" ;;
        real_pangenomes|simulated_pangenomes_perfect|simulated_pangenomes_flip)
                 dataset="$2" ;;
        *)
          echo "Invalid dataset: ${2:-<missing>}"
          echo "Allowed values: real, perfect, flip"
          exit 1
          ;;
      esac
      shift 2
      ;;
    *)
      echo "Usage: $0 --dataset <real|perfect|flip>"
      exit 1
      ;;
  esac
done

# Require dataset flag
if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real|perfect|flip>"
    exit 1
fi

goldfinder_dir="$(realpath "${dataset}/goldfinder_runs")"
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs")"

# Safety: bail if directory isn’t found
for dir in "$goldfinder_dir" "$coinfinder_dir"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in goldfinder_runs
for run_dir in "$goldfinder_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    nodes_file="$coinfinder_dir/$run_id/coincident_nodes_all.tsv"
    d_cutoff_file="$coinfinder_dir/$run_id/d_cutoff/${run_id}_d_cutoff.txt"
    pairs_file="$goldfinder_dir/$run_id/simultaneous_association_significant_pairs.csv"

    # Check that required files exist
    missing=false
    for f in "$nodes_file" "$d_cutoff_file" "$pairs_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Run d_distribution_goldfinder.R
    echo "  Running d_distribution_goldfinder.R for $run_id..."
    Rscript "$SCRIPT_DIR/d_distribution_goldfinder.R" \
        "$(realpath "$nodes_file")" \
        "$(realpath "$d_cutoff_file")" \
        "$(realpath "$pairs_file")"

    echo "  Finished $run_id"
    echo
done
