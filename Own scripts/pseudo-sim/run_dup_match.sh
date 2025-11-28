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
      echo "Usage: $0 --dataset <perfect|flip> [--mode <unfiltered|filtered>] [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <perfect|flip>"
    exit 1
fi

# Paths
if [ "$scope" = "all" ]; then
    gpa_dir="${dataset}/gpa_matches_all"
    summary_file="${dataset}/${dataset}_dup_match_${mode}_all.tsv"
    match_file="real_pangenomes/matched_all.csv"
else
    gpa_dir="${dataset}/gpa_matches"
    summary_file="${dataset}/${dataset}_dup_match_${mode}_selected.tsv"
    match_file="real_pangenomes/species_categories.csv"
fi

script="scripts/pseudo-sim/dup_match.py"

# Overwrite the summary file
: > "$summary_file"

run_tool () {
    local base_runs_subdir=$1
    local subfolder=$2
    local runs_subdir="${base_runs_subdir}_${scope}"
    local runs_dir="${dataset}/${runs_subdir}"

    echo "Processing ${runs_subdir} (${subfolder}/${base_runs_subdir})"

    shopt -s nullglob
    local infiles=()

    if [[ "$base_runs_subdir" == "panforest_runs" ]]; then
        runs_dir="${runs_dir}/${mode}"
        infiles=( "$runs_dir"/*/"${subfolder}"/*dvalues_*.csv )
    else
        infiles=( "$runs_dir"/*/"${subfolder}"/*dvalues_*.csv )
    fi
    shopt -u nullglob

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
            -m "$match_file" \
            -o "$summary_file"
    done
}

# Run tools
run_tool "coinfinder_runs" "d_cutoff"
run_tool "goldfinder_runs" "d_distribution"
run_tool "panforest_runs" "imp_cutoff"

echo "Summary written to $summary_file"
