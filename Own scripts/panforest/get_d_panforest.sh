#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments
dataset=""
mode="unfiltered"
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
    --mode)
      case ${2:-} in
        unfiltered|filtered) mode="$2" ;;
        *)
          echo "Invalid mode: ${2:-<missing>}"
          echo "Allowed values: unfiltered, filtered"
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
      echo "Usage: $0 --dataset <real|perfect|flip> [--mode <unfiltered|filtered>] [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real|perfect|flip>"
    exit 1
fi

# Scope-aware directories
panforest_dir="$(realpath "${dataset}/panforest_runs_${scope}/${mode}")"
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs_${scope}")"

# Choose GPA directory based on scope (not directly used here, but consistent)
if [ "$scope" = "all" ]; then
    gpa_dir="$(realpath "${dataset}/gpa_matches_all_not_pruned")"
else
    gpa_dir="$(realpath "${dataset}/gpa_matches")"
fi

# Safety: bail if directory isn’t found
for d in "$panforest_dir" "$coinfinder_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each run directory
for run_dir in "$panforest_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id (${mode}, scope=${scope})"

    imp_file="${run_dir}/imp.csv"
    perf_file="${run_dir}/performance.csv"
    dval_file="${run_dir}imp_cutoff/${run_id}_nodes.tsv"

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

    echo "  Running imp_distribution.R for $run_id..."
    Rscript "$SCRIPT_DIR/imp_distribution.R" \
        "$(realpath "$imp_file")" \
        "$(realpath "$perf_file")" \
        "$(realpath "$dval_file")" \
        "$dcutoff_value"

    echo "  Finished $run_id (${mode}, scope=${scope})"
    echo
done
