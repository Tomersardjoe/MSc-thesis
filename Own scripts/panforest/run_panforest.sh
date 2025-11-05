#!/usr/bin/env bash

required_env="panforest"
if [ "${CONDA_DEFAULT_ENV:-}" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Defaults
dataset=""
mode="unfiltered"
scope="selected"
MAX_JOBS=4

# Parse arguments
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

# Directories depend on scope
if [ "$scope" = "selected" ]; then
    gpa_dir="${dataset}/gpa_matches"
    out_base="${dataset}/panforest_runs_selected/${mode}"
    parallel=false
    panforest_args=(-n 1000 -d 16)
else
    gpa_dir="${dataset}/gpa_matches_all"
    out_base="${dataset}/panforest_runs_all/${mode}"
    parallel=true
    panforest_args=(-n 500 -d 8)
fi

pan_dir="panforest"
script_dir="scripts/panforest"

if [ ! -d "$gpa_dir" ]; then
    echo "Error: Directory '$gpa_dir' not found. Check the path."
    exit 1
fi

mkdir -p "$out_base"

# Preprocess matrices
for gpa_file in "$gpa_dir"/*_REDUCED.csv; do
    filename=$(basename "$gpa_file")
    base="${filename%_REDUCED.csv}"
    species_taxid="${base##*_}"

    echo "[$(date +%H:%M:%S)] Preprocessing ${species_taxid} (${mode}, scope=${scope})..."

    outdir="${out_base}/${species_taxid}"
    mkdir -p "$outdir"

    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    if [ "$mode" = "unfiltered" ]; then
        python3 "$pan_dir/process_matrix.py" \
            -i "$gpa_file" \
            -o "${outdir}/collapsed_matrix.csv" \
            -d "$outdir"
    else
        dstat_file="${dataset}/coinfinder_runs_${scope}/${species_taxid}/coincident_nodes_all.tsv"
        if [ ! -f "$dstat_file" ]; then
            echo "Warning: D-stat file not found for ${species_taxid}, skipping."
            continue
        fi
        filtered_gpa="${outdir}/${species_taxid}_filtered.csv"
        python3 "$script_dir/pre_forest_d_filter.py" \
            "$gpa_file" \
            "$dstat_file" \
            -o "$filtered_gpa"
        python3 "$pan_dir/process_matrix.py" \
            -i "$filtered_gpa" \
            -o "${outdir}/collapsed_matrix.csv" \
            -d "$outdir"
    fi
done

# Run PanForest
for matrix in "$out_base"/*/collapsed_matrix.csv; do
    species_taxid=$(basename "$(dirname "$matrix")")
    outdir=$(dirname "$matrix")

    run_cmd() {
      echo "[$(date +%H:%M:%S)] Starting PanForest for ${species_taxid} (${mode}, scope=${scope})..."
      NUMBA_NUM_THREADS=24 python3 "$pan_dir/PanForest.py" \
          "${panforest_args[@]}" \
          -m "$matrix" \
          -pres 1 \
          -abs 1 \
          -o "$outdir"
      echo "[$(date +%H:%M:%S)] Finished PanForest for ${species_taxid} (${mode}, scope=${scope})"
    }

    if [ "$parallel" = true ]; then
        run_cmd &
        if (( $(jobs -r | wc -l) >= MAX_JOBS )); then
            wait -n
        fi
    else
        run_cmd
    fi
done

if [ "$parallel" = true ]; then
    wait
    echo "All panforest runs completed."
fi
