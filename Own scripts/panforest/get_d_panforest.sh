#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

panforest_dir="$(realpath "${dataset}/panforest_runs")"
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs")"

# Safety: bail if directory isn’t found
if [ ! -d "$panforest_dir" ]; then
    echo "Error: Directory '$panforest_dir' not found. Check the path."
    exit 1
fi

# Loop through each run directory
for run_dir in "$panforest_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    imp_file="${run_dir}/imp.csv"
    perf_file="${run_dir}/performance.csv"
    dval_file="${run_dir}imp_cutoff/${run_id}_nodes.tsv"

    # Check that all required files exist
    missing=false
    for f in "$imp_file" "$perf_file" "$dval_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            echo
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    dcutoff_file="${coinfinder_dir}/${run_id}/d_cutoff/${run_id}_d_cutoff.txt"
    
    if [ -f "$dcutoff_file" ]; then
        dcutoff_value=$(<"$dcutoff_file")
    else
        echo "  No D-value cutoff found for $run_id, will be calculated in R script"
        dcutoff_value="NA"
    fi

    # Run imp_distribution.R
    echo "  Running imp_distribution.R for $run_id..."
    Rscript "$SCRIPT_DIR/imp_distribution.R" \
        "$(realpath "$imp_file")" \
        "$(realpath "$perf_file")" \
        "$(realpath "$dval_file")" \
        "$dcutoff_value"

    echo "  Finished $run_id"
    echo
done
