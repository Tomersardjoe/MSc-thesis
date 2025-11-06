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
scope="selected"   # default
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
      echo "Usage: $0 --dataset <real|perfect|flip> [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

# Require dataset flag
if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real|perfect|flip>"
    exit 1
fi

# Scope-aware coinfinder directory
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs_${scope}")"

# Choose GPA directory based on scope
if [ "$scope" = "all" ]; then
    gpa_dir="$(realpath "${dataset}/gpa_matches_all")"
else
    gpa_dir="$(realpath "${dataset}/gpa_matches")"
fi

# Safety: bail if directories aren’t found
for d in "$coinfinder_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in coinfinder_runs_${scope}
for run_dir in "$coinfinder_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id (scope=${scope})"

    nodes_file="${run_dir}coincident_nodes_all.tsv"
    pairs_file="${run_dir}coincident_pairs.tsv"
    tree_file="${run_dir}/${run_id}_fixed.nwk"
    gpa_file=$(ls "${gpa_dir}"/*"${run_id}"_REDUCED.csv 2>/dev/null | head -n1)

    # Check that all required files exist
    missing=false
    for f in "$nodes_file" "$pairs_file" "$tree_file" "$gpa_file"; do
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
        "$(realpath "$gpa_file")"

    echo "  Finished $run_id (scope=${scope})"
    echo
done
