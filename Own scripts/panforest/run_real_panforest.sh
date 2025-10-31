#!/usr/bin/env bash

required_env="panforest"
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
    *)
      echo "Usage: $0 --dataset <real|perfect|flip> [--mode <unfiltered|filtered>]"
      exit 1
      ;;
  esac
done

if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real|perfect|flip>"
    exit 1
fi

gpa_dir="${dataset}/gpa_matches"
pan_dir="panforest"
script_dir="scripts/panforest"
out_base="${dataset}/panforest_runs/${mode}"

if [ ! -d "$gpa_dir" ]; then
    echo "Error: Directory '$gpa_dir' not found. Check the path."
    exit 1
fi

mkdir -p "$out_base"

for gpa_file in "$gpa_dir"/*_REDUCED.csv; do
    filename=$(basename "$gpa_file")
    base="${filename%_REDUCED.csv}"
    species_taxid="${base##*_}"

    echo "Starting ${species_taxid} (${mode})..."

    outdir="${out_base}/${species_taxid}"
    mkdir -p "$outdir"

    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    if [ "$mode" = "unfiltered" ]; then
        # Step 1: Process matrix
        python3 "$pan_dir/process_matrix.py" \
            -i "$gpa_file" \
            -o "${outdir}/collapsed_matrix.csv" \
            -d "$outdir"

        matrix="${outdir}/collapsed_matrix.csv"
    else
        # Step 1a: Filter GPA matrix by D > 0
        dstat_file="${dataset}/coinfinder_runs/${species_taxid}/coincident_nodes_all.tsv"
        if [ ! -f "$dstat_file" ]; then
            echo "Warning: D-stat file not found for ${species_taxid}, skipping."
            continue
        fi
        filtered_gpa="${outdir}/${species_taxid}_filtered.csv"
        python3 "$script_dir/pre_forest_d_filter.py" \
            "$gpa_file" \
            "$dstat_file" \
            -o "$filtered_gpa"
    
        # Step 1b: Process the filtered GPA file
        python3 "$pan_dir/process_matrix.py" \
            -i "$filtered_gpa" \
            -o "${outdir}/collapsed_matrix.csv" \
            -d "$outdir"
    
        matrix="${outdir}/collapsed_matrix.csv"
    fi

    # Step 2: Run PanForest
    python3 "$pan_dir/PanForest.py" \
        -n 1000 \
        -d 16 \
        -m "$matrix" \
        -pres 1 \
        -abs 1 \
        -o "$outdir"

    echo "Finished ${species_taxid} (${mode})"
done
