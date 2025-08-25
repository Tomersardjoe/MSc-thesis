#!/usr/bin/env bash

PARENT_DIR="${1:-}"
SIG_THRESH="${2:-0.05}"

if [ -z "$PARENT_DIR" ]; then
  >&2 echo "Exception: No Coinfinder directory argument passed to the script. Please provide a valid path."
  exit 1
fi

if [ ! -d "$PARENT_DIR" ]; then
  >&2 echo "Exception: '$PARENT_DIR' is not a valid directory."
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for type in closed moderate open; do
  for i in {1..10}; do
    echo "Starting $type rep${i}..."

    dataset_dir="${PARENT_DIR}/${type}/rep${i}"

    if [ ! -d "$dataset_dir" ]; then
      echo "Warning: $dataset_dir does not exist. Skipping."
      continue
    fi

    # Run R script
    Rscript "$SCRIPT_DIR/p_adj_coinfinder.R" "$dataset_dir" "$SIG_THRESH"

    # Run Python script
    python3 "$SCRIPT_DIR/filter_sig.py" "$dataset_dir" "$SIG_THRESH"

    echo "Finished $type rep${i}"
    echo "-----------------------------"
  done
done

echo "FDR-BH applied to all Coinfinder datasets."

