#!/usr/bin/env bash

required_env="coinfinder_new"
if [ "${CONDA_DEFAULT_ENV:-}" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Parse arguments
dataset=""
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
      echo "Usage: $0 --dataset <real|perfect|flip> [--scope <selected|all>]"
      exit 1
      ;;
  esac
done

# Require dataset flag
if [ -z "$dataset" ]; then
    echo "Error: You must provide --dataset <real|perfect|flip>"
    exit 1
fi

# Directories depend on scope
if [ "$scope" = "selected" ]; then
    tree_dir="real_pangenomes/tree_matches"
    gpa_dir="${dataset}/gpa_matches"
    out_base="${dataset}/coinfinder_runs_selected"
    extra_args=()
else
    tree_dir="real_pangenomes/tree_matches_all"
    gpa_dir="${dataset}/gpa_matches_all"
    out_base="${dataset}/coinfinder_runs_all"
    extra_args=(-x 4)
fi

# Safety: bail if directories aren’t found
for d in "$tree_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

mkdir -p "$out_base"

for tree_file in "$tree_dir"/*.nwk; do
    filename=$(basename "$tree_file")
    basename_noext="${filename%.nwk}"
    species_taxid="${basename_noext#reduced_}"

    gpa_file=$(ls "${gpa_dir}"/*"${species_taxid}"_REDUCED.tab 2>/dev/null | head -n1)
    if [ ! -f "${gpa_file:-}" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    echo "Starting ${species_taxid} (scope=${scope})..."

    outdir="${out_base}/${species_taxid}"
    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "Skipping ${species_taxid}: $outdir already exists and is not empty."
        continue
    fi

    gpa_file="$(realpath "$gpa_file")"
    tree_file="$(realpath "$tree_file")"
    mkdir -p "$outdir"

    fixed_tree="${outdir}/${species_taxid}_fixed.nwk"
    Rscript -e "
      suppressMessages(library(ape));
      tr <- read.tree('$tree_file');
      tr\$edge.length[tr\$edge.length <= 1e-8] <- 1e-6;
      write.tree(tr, file='$fixed_tree')
    "

    (
      cd "$outdir" || exit
      NUMBA_NUM_THREADS=24 coinfinder \
        -i "$gpa_file" \
        -p "${species_taxid}_fixed.nwk" \
        -a \
        -n \
        -E \
        "${extra_args[@]}"
    )

    echo "Finished ${species_taxid}"
done
