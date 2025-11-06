#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

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
      mode="$2"
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

# Locate directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
ROOT_DIR="$(dirname "$PARENT_DIR")"

# Top-level output directory under repo root
base_outdir="${ROOT_DIR}/combined_results/${dataset}_${mode}_${scope}"
mkdir -p "$base_outdir"

# Dup summary file (scope-aware)
if [ "$dataset" = "real_pangenomes" ]; then
    summary_file="$(realpath "${dataset}/${dataset}_match_${mode}_${scope}.tsv")"
    if [ ! -f "$summary_file" ]; then
        echo "Error: ${dataset}_match_${mode}_${scope}.tsv not found. Exiting."
        exit 1
    fi
    sum_file="$summary_file"
else
    dup_summary_file="$(realpath "${dataset}/${dataset}_dup_match_${mode}_${scope}.tsv")"
    if [ ! -f "$dup_summary_file" ]; then
        echo "Error: dup_match_${mode}_${scope}.tsv not found for dataset $dataset. Exiting."
        exit 1
    fi
    sum_file="$dup_summary_file"
fi

# Run the R script
echo "Running confusion_plots_all.R for dataset '$dataset' (mode: $mode, scope: $scope)"
cmd=( Rscript "$SCRIPT_DIR/confusion_plots_all.R" \
      "$sum_file" \
      "$base_outdir" )

"${cmd[@]}"

echo "Finished generating confusion plots for dataset '$dataset' (mode: $mode, scope: $scope)"
