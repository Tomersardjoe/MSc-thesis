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

coinfinder_dir="$(realpath "${dataset}/coinfinder_runs_${scope}")"
goldfinder_dir="$(realpath "${dataset}/goldfinder_runs_${scope}")"
panforest_dir="$(realpath "${dataset}/panforest_runs_${scope}/${mode}")"

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
    echo "Processing run: $run_id (${mode}, scope=${scope})"

    coin_file="${coinfinder_dir}/${run_id}/d_cutoff/coinfinder_dvalues_${run_id}.csv"
    gold_file="${goldfinder_dir}/${run_id}/d_distribution/goldfinder_dvalues_${run_id}.csv"
    pan_file="${panforest_dir}/${run_id}/imp_cutoff/panforest_dvalues_${run_id}.csv"
    d_cutoff_file="${coinfinder_dir}/${run_id}/d_cutoff/${run_id}_d_cutoff.txt"
    
    # Check that required cutoff file exists
    if [ ! -f "$d_cutoff_file" ]; then
        echo "  Skipping $run_id - missing file: $d_cutoff_file"
        continue
    fi
    
    # Count how many files exist and have data
    has_data() {
        local f="$1"
        [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]
    }
    
    count=0
    has_data "$coin_file" && count=$((count+1))
    has_data "$gold_file" && count=$((count+1))
    has_data "$pan_file"  && count=$((count+1))
    
    # If all three missing, bail
    if [ $count -eq 0 ]; then
        echo "Error: all three files are missing for $run_id. Skipping."
        continue
    fi
    
    # In perfect dataset, if only PanForest is present, bail
    if [ "$dataset" = "simulated_pangenomes_perfect" ] && has_data "$pan_file" && \
       ! has_data "$coin_file" && ! has_data "$gold_file"; then
        echo "Error: perfect dataset run $run_id has only PanForest data. Exiting."
        continue
    fi
        
    # Run-specific output directory
    run_outdir="${base_outdir}/${run_id}"
    mkdir -p "$run_outdir"

    # Build the command
    cmd=( Rscript "$SCRIPT_DIR/comb_plots.R"
          "$(realpath "$coin_file")"
          "$(realpath "$gold_file")"
          "$(realpath "$pan_file")"
          "$(realpath "$d_cutoff_file")" )

    # If dataset is simulated, require the dup_summary file
    if [[ "$dataset" == simulated_pangenomes_* ]]; then
        dup_summary_file="$(realpath "${dataset}/${dataset}_dup_match_${mode}_${scope}.tsv")"
        if [ -f "$dup_summary_file" ]; then
            cmd+=( "$dup_summary_file" )
        else
            echo "Error:${dataset}/${dataset}_dup_match_${mode}_${scope}.tsv not found for dataset $dataset. Exiting."
            exit 1
        fi
    fi

    # Append run-specific output directory
    cmd+=( "$run_outdir" )

    # Run the R script with built command
    echo "  Running comb_plots.R for $run_id..."
    "${cmd[@]}"

    echo "  Finished $run_id (${mode}, scope=${scope})"
    echo
done
