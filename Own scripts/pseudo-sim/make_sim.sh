#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

dup_mode="perfect"   # default
scope="selected"     # default

while [[ $# -gt 0 ]]; do
  case $1 in
    --dup_mode)
      dup_mode="$2"
      shift 2
      ;;
    --dup_mode=*)
      dup_mode="${1#*=}"
      shift
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
    --scope=*)
      scope="${1#*=}"
      if [[ "$scope" != "selected" && "$scope" != "all" ]]; then
        echo "Invalid scope: $scope"
        echo "Allowed values: selected, all"
        exit 1
      fi
      shift
      ;;
    *)
      echo "Usage: $0 [--dup_mode <perfect|flip>] [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

# Scope-aware directories
nodes_dir="real_pangenomes/coinfinder_runs_${scope}"
gpa_dir="real_pangenomes/gpa_matches"
script="scripts/pseudo-sim/sim_pan.py"
out_base="simulated_pangenomes_${dup_mode}"

# Safety checks
for d in "$nodes_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found."
        exit 1
    fi
done

mkdir -p "$out_base"

# Loop through each nodes_all.tsv file
for nodes_file in "$nodes_dir"/*/coincident_nodes_all.tsv; do
    species_taxid=$(basename "$(dirname "$nodes_file")")

    gpa_file="$gpa_dir/${species_taxid}_REDUCED.csv"
    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid} (scope=${scope})..."

    outdir="${out_base}"
    mkdir -p "$outdir/gpa_matches"

    nodes_file="$(realpath "$nodes_file")"
    gpa_file="$(realpath "$gpa_file")"
    outdir="$(realpath "$outdir")"

    python3 "$script" \
        "$nodes_file" \
        "$gpa_file" \
        -o "${outdir}/gpa_matches/duplicates_${species_taxid}_REDUCED.csv" \
        --tab_out "${outdir}/gpa_matches/duplicates_${species_taxid}_REDUCED.tab" \
        --dup_mode "$dup_mode"

    echo "Finished ${species_taxid}"
done
