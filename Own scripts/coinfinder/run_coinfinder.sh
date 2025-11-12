#!/usr/bin/env bash

required_env="coinfinder"
if [ "${CONDA_DEFAULT_ENV:-}" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Defaults
dataset=""
scope="selected"
MAX_JOBS=2

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
    parallel=false
else
    tree_dir="real_pangenomes/tree_matches_all"
    gpa_dir="${dataset}/gpa_matches_all"
    out_base="${dataset}/coinfinder_runs_all"
    extra_args=()
    parallel=true
fi

for d in "$tree_dir" "$gpa_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

mkdir -p "$out_base"

# Preprocess trees
for tree_file in "$tree_dir"/*.nwk; do
    filename=$(basename "$tree_file")
    basename_noext="${filename%.nwk}"
    species_taxid="${basename_noext#reduced_}"

    outdir="${out_base}/${species_taxid}"
    mkdir -p "$outdir"

    orig_tree="$(realpath "$tree_file")"
    fixed_tree="${outdir}/${species_taxid}_fixed.nwk"

    echo "[$(date +%H:%M:%S)] Preprocessing tree for ${species_taxid}..."
    Rscript -e "
      suppressMessages(library(ape));
      tr <- read.tree('$orig_tree');
      tr\$edge.length[tr\$edge.length <= 1e-8] <- 1e-6;
      write.tree(tr, file='$fixed_tree')
    "
done

# Run Coinfinder
count=0
for fixed_tree in "$out_base"/*/*_fixed.nwk; do
    species_taxid=$(basename "$fixed_tree" _fixed.nwk)
    outdir=$(dirname "$fixed_tree")

    gpa_file=$(ls "${gpa_dir}"/*"${species_taxid}"_REDUCED.tab 2>/dev/null | head -n1)
    if [ ! -f "$gpa_file" ]; then
        echo "Warning: No GPA file for ${species_taxid}, skipping."
        continue
    fi

    if [ -f "${outdir}/coincident_heatmap0.pdf" ]; then
        echo "Skipping ${species_taxid}: run already completed."
        continue
    fi

    gpa_file="$(realpath "$gpa_file")"
    fixed_tree_abs="$(realpath "$fixed_tree")"

    run_cmd() {
      echo "[$(date +%H:%M:%S)] Starting Coinfinder for ${species_taxid} (scope=${scope})..."
      (
        cd "$outdir" || exit
        nice -n 10 coinfinder \
          -i "$gpa_file" \
          -p "$fixed_tree_abs" \
          -a \
          -x 1 \
          -n \
          -E \
          "${extra_args[@]}"
        rm -f *.gexf *.gml # Remove network files
        echo "[$(date +%H:%M:%S)] Finished Coinfinder for ${species_taxid}"
      )
    }

    if [ "$parallel" = true ]; then
        run_cmd &
        ((++count % MAX_JOBS == 0)) && wait
    else
        run_cmd
    fi

done

if [ "$parallel" = true ]; then
    wait
    echo "All Coinfinder runs completed."
fi
