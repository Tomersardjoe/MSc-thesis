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
base_outdir="${ROOT_DIR}/combined_results/${dataset}_${mode}_${scope}/results"
mkdir -p "$base_outdir"

# Tool run directories (scope-aware)
coin_dir="$(realpath "${dataset}/coinfinder_runs${scope:+_$scope}")"
gold_dir="$(realpath "${dataset}/goldfinder_runs${scope:+_$scope}")"
pan_dir="$(realpath "${dataset}/panforest_runs${scope:+_$scope}/${mode}")"

for d in "$coin_dir" "$gold_dir" "$pan_dir"; do
  if [ ! -d "$d" ]; then
    echo "Error: Expected directory $d not found. Exiting."
    exit 1
  fi
done

# Category file
categories_file="$(realpath "real_pangenomes/species_categories.csv")"
if [ ! -f "$categories_file" ]; then
  echo "Error: species_categories.csv not found in $dataset. Exiting."
  exit 1
fi

# Run the R script
echo "Running aggregated_upset.R for dataset '$dataset' (mode: $mode, scope: $scope)"
cmd=( Rscript "$SCRIPT_DIR/aggregated_upset.R" \
      "$coin_dir" \
      "$gold_dir" \
      "$pan_dir" \
      "$categories_file" \
      "$base_outdir" )

"${cmd[@]}"

echo "Finished generating aggregated UpSet plots for dataset '$dataset' (mode: $mode, scope: $scope)"
echo "Output files created in $base_outdir"
