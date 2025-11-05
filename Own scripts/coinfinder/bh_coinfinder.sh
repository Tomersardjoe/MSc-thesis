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

# Loop over all subdirectories
for dataset_dir in "$PARENT_DIR"/*/; do
  # Skip if not a directory
  [ -d "$dataset_dir" ] || continue

  echo "Starting $(basename "$dataset_dir")..."

  # Run R script
  Rscript "$SCRIPT_DIR/p_adj_coinfinder.R" "$dataset_dir" "$SIG_THRESH"

  # Run Python script
  python3 "$SCRIPT_DIR/filter_sig.py" "$dataset_dir" "$SIG_THRESH"

  echo "Finished $(basename "$dataset_dir")"
  echo "-----------------------------"
done

echo "FDR-BH applied to all Coinfinder subdirectories."
