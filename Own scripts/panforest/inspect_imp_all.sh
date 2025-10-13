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

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

panforest_dir="$(realpath "${dataset}/panforest_runs")"
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs")"
panscript_dir="$(realpath "panforest")"

# Safety: bail if directories aren’t found
for d in "$panforest_dir" "$coinfinder_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in panforest_runs
for run_dir in "$panforest_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    imp_file="${run_dir}/imp.csv"
    nodes_file="${run_dir}imp_cutoff/${run_id}_nodes.tsv"
    cutoff_file="${run_dir}imp_cutoff/${run_id}_cutoff_value.txt"
    dcutoff_file="${coinfinder_dir}/${run_id}/d_cutoff/${run_id}_d_cutoff.txt"

    # Check that all required files exist
    missing=false
    for f in "$imp_file" "$nodes_file" "$cutoff_file" "$dcutoff_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Read cutoff values from files
    cutoff_value=$(<"$cutoff_file")
    dcutoff_value=$(<"$dcutoff_file")
    
    # Simplify importance matrix
    echo "  Running simplify_imp.py for $run_id..."
    python3 "$panscript_dir/simplify_imp.py" \
        "$cutoff_value" \
        "$(realpath "$imp_file")" \
        "${run_dir}imp_cutoff/imp_simplified.csv"
    
    # Run inspect_imp.R
    echo "  Running inspect_imp.R for $run_id..."
    Rscript "$SCRIPT_DIR/inspect_imp.R" \
        "$(realpath "$imp_file")" \
        "$(realpath "${run_dir}imp_cutoff/imp_simplified.csv")" \
        "$(realpath "$nodes_file")" \
        "$cutoff_value" \
        "$dcutoff_value"

    echo "  Finished $run_id"
    echo
done
