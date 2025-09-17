#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

coinfinder_dir="$(realpath "real_pangenomes/coinfinder_runs")"
tree_dir="$(realpath "real_pangenomes/tree_matches")"
pa_dir="$(realpath "real_pangenomes/gpa_matches")"

# Safety: bail if directories aren’t found
for d in "$coinfinder_dir" "$tree_dir" "$pa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in coinfinder_runs
for run_dir in "$coinfinder_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"

    nodes_file="${run_dir}coincident_nodes.tsv"
    pairs_file="${run_dir}coincident_pairs.tsv"
    tree_file="${tree_dir}/reduced_${run_id}.nwk"
    pa_file="${pa_dir}/${run_id}_REDUCED.csv"

    # Check that all required files exist
    missing=false
    for f in "$nodes_file" "$pairs_file" "$tree_file" "$pa_file"; do
        if [ ! -f "$f" ]; then
            echo "  Skipping $run_id - missing file: $f"
            missing=true
        fi
    done
    if [ "$missing" = true ]; then
        continue
    fi

    # Run d_distribution_coinfinder.R
    echo "  Running d_distribution_coinfinder.R for $run_id..."
    Rscript "$SCRIPT_DIR/d_distribution_coinfinder.R" \
        "$(realpath "$nodes_file")" \
        "$(realpath "$pairs_file")" \
        "$(realpath "$tree_file")" \
        "$(realpath "$pa_file")"

    echo "  Finished $run_id"
    echo
done
