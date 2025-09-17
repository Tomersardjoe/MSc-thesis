#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

coinfinder_dir="$(realpath "real_pangenomes/coinfinder_runs")"
goldfinder_dir="$(realpath "real_pangenomes/goldfinder_runs")"
panforest_dir="$(realpath "real_pangenomes/panforest_runs")"

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
    pan_file="${panforest_dir}/${run_id}/panforest_dvalues_${run_id}.csv"

    # Check that all required files exist
    missing=false
    for f in "$coin_file" "$gold_file" "$pan_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Run the comparison R script
    echo "  Running comp_plots.R for $run_id..."
    Rscript "$SCRIPT_DIR/comp_plots.R" \
        "$(realpath "$coin_file")" \
        "$(realpath "$gold_file")" \
        "$(realpath "$pan_file")"

    echo "  Finished $run_id"
    echo
done
