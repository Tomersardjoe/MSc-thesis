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
      if [[ "${2:-}" == "real" ]]; then
        dataset="real_pangenomes"
      else
        echo "Invalid dataset: ${2:-<missing>}"
        echo "Allowed value: real"
        exit 1
      fi
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
      echo "Usage: $0 --dataset real [--mode <unfiltered|filtered>] [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset real"
    exit 1
fi

# Paths
if [ "$scope" = "all" ]; then
    gpa_dir="${dataset}/gpa_matches_all"
    match_file="${dataset}/matched_all.csv"
else
    gpa_dir="${dataset}/gpa_matches"
    match_file="${dataset}/matched_selected.csv"
fi

script="scripts/combined/real_match.py"
summary_file="${dataset}/${dataset}_match_${mode}_${scope}.tsv"

# Safety checks
if [ ! -d "$gpa_dir" ]; then
    echo "Error: GPA directory '$gpa_dir' not found."
    exit 1
fi

if [ ! -f "$match_file" ]; then
    echo "Error: Match file '$match_file' not found."
    exit 1
fi

# Overwrite the summary file
: > "$summary_file"

run_tool () {
    local base_runs_subdir=$1
    local subfolder=$2
    local runs_subdir="${base_runs_subdir}_${scope}"
    local runs_dir="${dataset}/${runs_subdir}"

    # If panforest, include mode
    if [[ "$base_runs_subdir" == "panforest_runs" ]]; then
        runs_dir="${runs_dir}/${mode}"
    fi

    echo "Processing ${runs_subdir} (${subfolder}${base_runs_subdir=="panforest_runs" ? " | mode: ${mode}" : ""})"

    shopt -s nullglob
    local infiles=( "$runs_dir"/*/${subfolder}/*dvalues_*.csv )
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

        echo "Run ${count} / ${total} (pangenome ${run_id})"
        python3 "$script" \
            -i "$infile" \
            -a "$match_file" \
            -o "$summary_file"
    done
}

# Run tools
found_any=false

[ -d "${dataset}/coinfinder_runs_${scope}" ] && { run_tool "coinfinder_runs" "d_cutoff"; found_any=true; }
[ -d "${dataset}/goldfinder_runs_${scope}" ] && { run_tool "goldfinder_runs" "d_distribution"; found_any=true; }
[ -d "${dataset}/panforest_runs_${scope}" ] && { run_tool "panforest_runs" "imp_cutoff"; found_any=true; }

if [ "$found_any" = false ]; then
    echo "Error: No run directories found under $dataset."
    exit 1
fi

echo "Summary written to $summary_file"
