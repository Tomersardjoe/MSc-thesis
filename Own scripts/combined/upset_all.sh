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
if [ "$scope" = "selected" ]; then
  categories_file="$(realpath "real_pangenomes/species_categories.csv")"
elif [ "$scope" = "all" ]; then
  categories_file="$(realpath "real_pangenomes/species_categories.csv")"
else
  echo "Error: Unknown scope '$scope'"
  exit 1
fi

if [ ! -f "$categories_file" ]; then
  echo "Error: species_categories.csv not found at $categories_file. Exiting."
  exit 1
fi

# Build the base command
cmd=( Rscript "$SCRIPT_DIR/aggregated_upset.R" \
      "$coin_dir" \
      "$gold_dir" \
      "$pan_dir" \
      "$categories_file" \
      "$base_outdir" )

# Optional summary file if dataset is simulated
if [[ "$dataset" == "simulated_pangenomes_flip" || "$dataset" == "simulated_pangenomes_perfect" ]]; then
  dataset_suffix="${dataset##simulated_pangenomes_}"
  summary_file="simulated_pangenomes_${dataset_suffix}_dup_match_${mode}_${scope}.tsv"
  summary_path="$(realpath "$dataset/$summary_file")"

  if [ -f "$summary_path" ]; then
    cmd+=( "$summary_path" )
    echo "Passed $summary_path as optional argument"
  else
    echo "Warning: expected summary file $summary_path not found, continuing without it."
  fi
fi

# Run the R script
echo "Running aggregated_upset.R for dataset '$dataset' (mode: $mode, scope: $scope)"
"${cmd[@]}"

echo "Finished generating aggregated UpSet plots for dataset '$dataset' (mode: $mode, scope: $scope)"
echo "Output files created in $base_outdir"
