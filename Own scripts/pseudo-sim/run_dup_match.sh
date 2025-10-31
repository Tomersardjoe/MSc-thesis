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

gpa_dir="${dataset}/gpa_matches"
script="scripts/pseudo-sim/dup_match.py"
summary_file="${dataset}/dup_match_${mode}.tsv"
species_categories_file="real_pangenomes/species_categories.csv"

# Overwrite the summary file
: > "$summary_file"

run_tool () {
    local runs_subdir=$1
    local subfolder=$2
    local runs_dir="${dataset}/${runs_subdir}"

    # If panforest, include mode
    if [[ "$runs_subdir" == "panforest_runs" ]]; then
        runs_dir="${runs_dir}/${mode}"
    fi

    # Safety checks
    for d in "$runs_dir" "$gpa_dir"; do
        if [ ! -d "$d" ]; then
            echo "Error: Directory '$d' not found."
            return
        fi
    done
    if [ ! -f "$species_categories_file" ]; then
        echo "Error: Species categories file '$species_categories_file' not found."
        return
    fi

    echo "Processing ${runs_subdir} (${subfolder}/${runs_subdir=="panforest_runs" ? " | mode: ${mode}" : ""})"

    shopt -s nullglob # if no match is found, empty string
    local infiles=( "$runs_dir"/*/${subfolder}/*dvalues_*.csv )
    shopt -u nullglob # turn off file matcher

    local total=${#infiles[@]}
    if (( total == 0 )); then
        echo "No input files found under ${runs_dir}/${subfolder}. Skipping."
        return
    fi

    local count=0
    for infile in "${infiles[@]}"; do
        count=$((count+1))
        local run_id
        run_id=$(basename "$infile" | grep -oE '[0-9]+')

        local dup_file="$gpa_dir/duplicates_${run_id}_REDUCED.csv"
        if [ ! -f "$dup_file" ]; then
            echo "Warning: No GPA file for run ${run_id}, skipping."
            continue
        fi

        echo "Run ${count} / ${total} (pangenome ${run_id})"
        python3 "$script" \
            -i "$infile" \
            -d "$dup_file" \
            -s "$species_categories_file" \
            -o "$summary_file"
    done
}

# Run tools
run_tool "coinfinder_runs" "d_cutoff"
run_tool "goldfinder_runs" "d_distribution"
run_tool "panforest_runs" "imp_cutoff"

echo "Summary written to $summary_file"
