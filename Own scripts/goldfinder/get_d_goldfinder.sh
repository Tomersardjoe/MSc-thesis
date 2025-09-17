#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

base_dir="real_pangenomes"
goldfinder_dir="$(realpath "$base_dir/goldfinder_runs")"
coinfinder_dir="$(realpath "$base_dir/coinfinder_runs")"

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
    pairs_file="$goldfinder_dir/$run_id/simultaneous_association_significant_pairs.csv"

    # Check that required files exist
    missing=false
    for f in "$nodes_file" "$pairs_file"; do
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
        "$(realpath "$pairs_file")"

    echo "  Finished $run_id"
    echo
done
