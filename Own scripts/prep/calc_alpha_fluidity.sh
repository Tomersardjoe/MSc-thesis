#!/bin/bash

required_env="analysis"
if [ "$CONDA_DEFAULT_ENV" != "$required_env" ]; then
    echo "Please activate the '$required_env' environment before running this script."
    exit 1
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --gpa)
      gpa_all_dir="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    *)
      echo "Usage: $0 --gpa <path/to/gpa_matches_all> --outdir <output_dir>"
      exit 1
      ;;
  esac
done

# Safety checks
if [[ -z "$gpa_all_dir" || -z "$outdir" ]]; then
  echo "Error: both --gpa and --outdir must be provided"
  exit 1
fi
if [[ ! -d "$gpa_all_dir" ]]; then
  echo "Error: GPA directory '$gpa_all_dir' not found"
  exit 1
fi
mkdir -p "$outdir"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_csv="${outdir}/alpha_fluidity_all.csv"

# Write header if not present
if [ ! -f "$output_csv" ]; then
  echo "species_taxid,filename,fluidity_calc,openness" > "$output_csv"
fi

# Process each *_REDUCED.csv file
for file in "$gpa_all_dir"/*_REDUCED.csv; do
  filename=$(basename "$file")
  species_taxid=$(echo "$filename" | cut -d'_' -f1)

  # Skip if already processed
  pattern=$(printf '%s' "$species_taxid" | sed -e 's/[]\/$*.^|[]/\\&/g')
  if sed 's/\r$//' "$output_csv" | grep -Eq "^$pattern,"; then
    echo "Skipping $filename..."
    continue
  fi

  echo "Processing $filename..."
  result=$(Rscript "$SCRIPT_DIR/fluidity_calc.R" "$file")

  fluidity=$(echo "$result" | grep "Genomic fluidity" | sed -E 's/.*: *([0-9.]+) *$/\1/')
  openness=$(echo "$result" | grep "Pangenome openness" | sed -E 's/.*: *([0-9.]+) *$/\1/')

  echo "$species_taxid,$filename,$fluidity,$openness" >> "$output_csv"
done

echo "All results written to $output_csv"

