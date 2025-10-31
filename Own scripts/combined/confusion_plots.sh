#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Parse arguments
dataset=""
mode="unfiltered"
while [[ $# -gt 0 ]]; do
  case $1 in
    --dataset)
      case ${2:-} in
        perfect) dataset="simulated_pangenomes_perfect" ;;
        flip)    dataset="simulated_pangenomes_flip" ;;
        simulated_pangenomes_perfect|simulated_pangenomes_flip)
                 dataset="$2" ;;
        *)
          echo "Invalid dataset: ${2:-<missing>}"
          echo "Allowed values: perfect, flip"
          exit 1
          ;;
      esac
      shift 2
      ;;
    --mode)
      mode="$2"
      shift 2
      ;;
    *)
      echo "Usage: $0 --dataset <perfect|flip> [--mode <unfiltered|filtered>]"
      exit 1
      ;;
  esac
done

if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <perfect|flip>"
    exit 1
fi

# Locate directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# Top-level output directory under repo root
base_outdir="${PARENT_DIR}/combined_results/${dataset}_${mode}"
mkdir -p "$base_outdir"

# Dup summary file
dup_summary_file="$(realpath "${dataset}/dup_match_${mode}.tsv")"
  if [ ! -f "$dup_summary_file" ]; then
    echo "Error: dup_match_${mode}.tsv not found for dataset $dataset. Exiting."
    exit 1
  fi

# Run the R script
echo "Running confusion_plots_all.R for dataset '$dataset' (mode: $mode)"
cmd=( Rscript "$SCRIPT_DIR/confusion_plots_all.R"
      "$dup_summary_file"
      "$base_outdir" )

"${cmd[@]}"

echo "Finished generating confusion plots for dataset '$dataset' (mode: $mode)"