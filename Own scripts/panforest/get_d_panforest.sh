#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
panforest_root="$(realpath "real_pangenomes/panforest_runs")"

# Safety: bail if directory isn’t found
if [ ! -d "$panforest_root" ]; then
    echo "Error: Directory '$panforest_root' not found. Check the path."
    exit 1
fi

# Loop through each run directory
for run_dir in "$panforest_root"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    imp_file="${run_dir}imp_cutoff/imp_fixed.csv"
    dval_file="${run_dir}imp_cutoff/${run_id}_nodes.tsv"

    # Check that both required files exist
    missing=false
    for f in "$imp_file" "$dval_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Run imp_distribution.R
    echo "  Running imp_distribution.R for $run_id..."
    Rscript "$SCRIPT_DIR/imp_distribution.R" \
        "$(realpath "$imp_file")" \
        "$(realpath "$dval_file")"

    echo "  Finished $run_id"
    echo
done
