#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Resolve script and project directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR"/.. && pwd)"
SUMMARIES_DIR="$PROJECT_ROOT/summaries"

# Create summaries directory if it doesn't exist already
mkdir -p "$SUMMARIES_DIR"

echo "[🔄] Combining tool summaries..."

# Define master header columns
MASTER_COLUMNS=(
  Tool
  Dataset
  Total_Pangenome_Genes
  Network_Nodes
  Network_Edges
  Possible_Gene_Pairs
  Association_Rate
  Avg_Degree
  Module_Count
  Avg_Genes_per_Module
  Modularity
)

# Write header to combined CSV
(
  IFS=,
  echo "${MASTER_COLUMNS[*]}"
) > "$SUMMARIES_DIR/combined_summary.csv"

# Function to merge a tool-specific summary into the master CSV
combine_by_header() {
  local TOOL="$1"
  local FILE="$SUMMARIES_DIR/$2"

  # Read header row into an array
  IFS=, read -r -a COLS < <(head -n1 "$FILE")

  # Build a map: column_name → index
  declare -A header_map
  for idx in "${!COLS[@]}"; do
    header_map["${COLS[idx]}"]=$idx
  done

  # Process each data row
  tail -n +2 "$FILE" | while IFS=, read -r -a vals; do
    # Map column values by header name
    declare -A row
    for col in "${COLS[@]}"; do
      row["$col"]="${vals[${header_map[$col]}]:-NA}"
    done

    # Assemble output fields in MASTER_COLUMNS order
    fields=(
      "$TOOL"
      "${row[Dataset]}"
      "${row[Total_Pangenome_Genes]}"
      "${row[Total_Networked_Genes]:-NA}"
      "${row[Significant_Associations]:-NA}"
      "${row[Possible_Gene_Pairs]:-NA}"
      "${row[Association_Rate]:-NA}"
      "${row[Avg_Degree]:-NA}"
      "${row[Module_Count]:-NA}"
      "${row[Avg_Genes_per_Module]:-NA}"
      "${row[Modularity]:-NA}"
    )

    # Append to combined_summary.csv
    ( IFS=,; echo "${fields[*]}" ) >> "$SUMMARIES_DIR/combined_summary.csv"
  done
}

# Merge both tool summaries
combine_by_header "Coinfinder" "coinfinder_summary.csv"
combine_by_header "Goldfinder"  "goldfinder_summary.csv"

echo "[✅] Created combined_summary.csv at $SUMMARIES_DIR/combined_summary.csv"
