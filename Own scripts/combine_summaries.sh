#!/usr/bin/env bash

echo "[🔄] Combining tool summaries..."

MASTER_COLUMNS=(\
  Tool Dataset Total_Pangenome_Genes Network_Nodes \
  Network_Edges Possible_Gene_Pairs Association_Rate \
  Avg_Degree Module_Count Avg_Genes_per_Module Modularity\
)

# write header
( IFS=,; echo "${MASTER_COLUMNS[*]}" ) > combined_summary.csv

combine_by_header() {
  local TOOL="$1"
  local FILE="$2"
  IFS=, read -r -a COLS < <(head -n1 "$FILE")

  # map column name → index
  declare -A header_map
  for idx in "${!COLS[@]}"; do
    header_map["${COLS[idx]}"]=$idx
  done

  tail -n +2 "$FILE" | while IFS=, read -r -a vals; do
    # build row map
    declare -A row
    for col in "${COLS[@]}"; do
      row["$col"]="${vals[${header_map[$col]}]}"
    done

    # assemble and print 1 line
    fields=(\
      "$TOOL" \
      "${row[Dataset]:-NA}" \
      "${row[Total_Pangenome_Genes]:-NA}" \
      "${row[Total_Networked_Genes]:-NA}" \
      "${row[Significant_Associations]:-NA}" \
      "${row[Possible_Gene_Pairs]:-NA}" \
      "${row[Association_Rate]:-NA}" \
      "${row[Avg_Degree]:-NA}" \
      "${row[Module_Count]:-NA}" \
      "${row[Avg_Genes_per_Module]:-NA}" \
      "${row[Modularity]:-NA}"\
    )

    ( IFS=,; echo "${fields[*]}" ) >> combined_summary.csv
  done
}

combine_by_header "Coinfinder" coinfinder_summary.csv
combine_by_header "Goldfinder"   goldfinder_summary.csv

echo "[✅] Combined_summary.csv created"
