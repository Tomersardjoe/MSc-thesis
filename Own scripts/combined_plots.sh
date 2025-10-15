#!/usr/bin/env bash

required_env="analysis"
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

# Locate directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

coinfinder_dir="$(realpath "${dataset}/coinfinder_runs")"
goldfinder_dir="$(realpath "${dataset}/goldfinder_runs")"
panforest_dir="$(realpath "${dataset}/panforest_runs")"

# Safety: bail if directories aren’t found
for d in "$coinfinder_dir" "$goldfinder_dir" "$panforest_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each run_id in coinfinder_runs
for run_dir in "$coinfinder_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    coin_file="${coinfinder_dir}/${run_id}/d_cutoff/coinfinder_dvalues_${run_id}.csv"
    gold_file="${goldfinder_dir}/${run_id}/d_distribution/goldfinder_dvalues_${run_id}.csv"
    pan_file="${panforest_dir}/${run_id}/imp_cutoff/panforest_dvalues_${run_id}.csv"
    d_cutoff_file="${coinfinder_dir}/${run_id}/d_cutoff/${run_id}_d_cutoff.txt"

    # Check that all required files exist
    missing=false
    for f in "$coin_file" "$gold_file" "$pan_file" "$d_cutoff_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Run the comparison R script with all arguments
    echo "  Running comp_plots.R for $run_id..."
    Rscript "$SCRIPT_DIR/comp_plots.R" \
        "$(realpath "$coin_file")" \
        "$(realpath "$gold_file")" \
        "$(realpath "$pan_file")" \
        "$(realpath "$d_cutoff_file")" \

    echo "  Finished $run_id"
    echo
done
