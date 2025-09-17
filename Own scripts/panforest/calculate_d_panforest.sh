#!/usr/bin/env bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

panforest_dir="$(realpath "real_pangenomes/panforest_runs")"
coinfinder_dir="$(realpath "real_pangenomes/coinfinder_runs")"

# Safety: bail if directories aren’t found
for d in "$panforest_dir" "$coinfinder_dir"; do
    if [ ! -d "$d" ]; then
        echo "Error: Directory '$d' not found. Check the path."
        exit 1
    fi
done

# Loop through each subdirectory in panforest_runs and fix family group labels
for run_dir in "$panforest_dir"/*/; do
    run_id=$(basename "$run_dir")
    echo "Processing run: $run_id"
    
    imp_csv="${run_dir}imp.csv"
    imp_fixed="${run_dir}imp_fixed.csv"

    # Check if imp.csv exists and imp_fixed.csv doesn't
    if [[ -f "$imp_csv" && ! -f "$imp_fixed" ]]; then
        sed 's/^\([^,]*\)_family_group\(,\|$\)/\1\2/' "$imp_csv" > "$imp_fixed"
        echo "  Fixed file created: $imp_fixed"
    elif [[ -f "$imp_fixed" ]]; then
        echo "  Skipping - $imp_fixed already exists."
    else
        echo "  No imp.csv found in $run_dir - skipping."
    fi

    cleaned_matrix="${run_dir}${run_id}_collapsed_matrix_clean.csv"

    # Skip if a _nodes.tsv already exists
    if ls "${run_dir}"*_nodes.tsv >/dev/null 2>&1; then
        echo "  Skipping $run_id - _nodes.tsv already exists."
        continue
    fi

    # Skip if cleaned matrix already exists
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
    Rscript "$SCRIPT_DIR/prep_d_calc.R" "$collapsed_matrix_abs"

    if [ ! -f "$cleaned_matrix" ]; then
        echo "  Error: Cleaned matrix not found for $run_id after prep_d_calc.R"
        continue
    fi

    # Step 2: Find matching fixed tree in coinfinder_runs
    fixed_tree="${coinfinder_dir}/${run_id}/${run_id}_fixed.nwk"
    if [ ! -f "$fixed_tree" ]; then
        echo "  Skipping $run_id - fixed tree not found: $fixed_tree"
        continue
    fi

    # Step 3: Run calculate_d.R
    (
      cd "$run_dir" || exit
      echo "  Running calculate_d.R for $run_id..."
      NUMBA_NUM_THREADS=24 Rscript "$SCRIPT_DIR/../panforest/calculate_d.R" \
          -a . \
          -t "$(realpath "$fixed_tree")" \
          -g "${run_id}_collapsed_matrix_clean.csv" \
          -c 1 \
          -o "$run_id"
    )

    echo "  Finished $run_id"
    echo
done
