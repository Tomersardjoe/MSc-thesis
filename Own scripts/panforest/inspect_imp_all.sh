#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

panforest_dir="$(realpath "real_pangenomes/panforest_runs")"
coinfinder_dir="$(realpath "real_pangenomes/coinfinder_runs")"

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

    imp_file="${run_dir}imp_fixed.csv"
    nodes_file="${run_dir}${run_id}_nodes.tsv"
    cutoff_file="${run_dir}${run_id}_cutoff_value.txt"
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

    # Run inspect_imp.R
    echo "  Running inspect_imp.R for $run_id..."
    Rscript "$SCRIPT_DIR/inspect_imp.R" \
        "$(realpath "$imp_file")" \
        "$(realpath "$nodes_file")" \
        "$(realpath "$cutoff_file")" \
        "$(realpath "$dcutoff_file")"

    echo "  Finished $run_id"
    echo
done
