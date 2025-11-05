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

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
panscript_dir="$(realpath "panforest")"
panforest_dir="$(realpath "${dataset}/panforest_runs_${scope}/${mode}")"
coinfinder_dir="$(realpath "${dataset}/coinfinder_runs_${scope}")"

# Safety: bail if directories aren’t found
for d in "$panforest_dir" "$coinfinder_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in panforest_runs_${scope}/<mode>
for run_dir in "$panforest_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id (${mode}, scope=${scope})"

    mkdir -p "${run_dir}imp_cutoff"

    imp_csv="${run_dir}imp.csv"
    cleaned_matrix="${run_dir}imp_cutoff/${run_id}_collapsed_matrix_clean.csv"
    nodes_tsv="${run_dir}imp_cutoff/${run_id}_nodes.tsv"

    if [ -f "$nodes_tsv" ]; then
        echo "  Skipping $run_id - $nodes_tsv already exists."
        continue
    fi

    if [ -f "$cleaned_matrix" ]; then
        echo "  Skipping $run_id - cleaned matrix already exists."
        continue
    fi

    collapsed_matrix="${run_dir}collapsed_matrix.csv"
    if [ ! -f "$collapsed_matrix" ]; then
        echo "  Skipping $run_id - no collapsed_matrix.csv found."
        continue
    fi

    # Step 1: Run prep_d_calc.R
    collapsed_matrix_abs="$(realpath "$collapsed_matrix")"
    echo "  Running prep_d_calc.R for $run_id..."
    Rscript "$SCRIPT_DIR/prep_d_calc.R" "$collapsed_matrix_abs" "${run_dir}imp_cutoff"

    if [ ! -f "$cleaned_matrix" ]; then
        echo "  Error: Cleaned matrix not found for $run_id after prep_d_calc.R"
        continue
    fi

    # Step 2: Find matching fixed tree in coinfinder_runs_${scope}
    fixed_tree="${coinfinder_dir}/${run_id}/${run_id}_fixed.nwk"
    if [ ! -f "$fixed_tree" ]; then
        echo "  Skipping $run_id - fixed tree not found: $fixed_tree"
        continue
    fi

    # Step 3: Run calculate_d*.R depending on scope
    if [ "$scope" = "selected" ]; then
        calc_script="calculate_d.R"
    else
        calc_script="calculate_d_100.R"
    fi

    (
      cd "$run_dir" || exit
      echo "  Running $calc_script for $run_id..."
      NUMBA_NUM_THREADS=24 Rscript "$panscript_dir/$calc_script" \
          -a . \
          -t "$(realpath "$fixed_tree")" \
          -g "imp_cutoff/${run_id}_collapsed_matrix_clean.csv" \
          -c 1 \
          -o "imp_cutoff"
    )

    echo "  Finished $run_id (${mode}, scope=${scope})"
    echo
done
