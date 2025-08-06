### 🧩 Combine Summaries with Unified Headers

echo "[🔄] Combining tool summaries..."

MASTER_COLUMNS=("Tool" "Dataset" "Total_Pangenome_Genes" "Network_Nodes" "Network_Edges" \
"Gene_Associations" "Possible_Gene_Pairs" "Association_Rate" "Avg_Degree" \
"Module_Count" "Avg_Genes_per_Module" "Modularity")

echo "${MASTER_COLUMNS[*]}" | tr ' ' ',' > combined_summary.csv

function combine_by_header() {
  TOOL=$1
  FILE=$2

  HEADER=$(head -n 1 "$FILE")
  IFS=',' read -ra COLS <<< "$HEADER"

  declare -A header_map
  for i in "${!COLS[@]}"; do
    header_map["${COLS[$i]}"]=$i
  done

  tail -n +2 "$FILE" | while IFS=',' read -ra values; do
    declare -A row_map
    for col in "${COLS[@]}"; do
      idx=${header_map["$col"]}
      row_map["$col"]="${values[$idx]}"
    done

    LINE="$TOOL"

    case "$TOOL" in
      "Coinfinder")
        LINE="$LINE,${row_map["Dataset"]:-NA}"
        LINE="$LINE,${row_map["Total_Pangenome_Genes"]:-NA}"
        LINE="$LINE,${row_map["Total_Networked_Genes"]:-NA}"    # Network_Nodes
        LINE="$LINE,NA"                                          # Network_Edges
        LINE="$LINE,${row_map["Significant_Associations"]:-NA}" # Gene_Associations
        LINE="$LINE,${row_map["Possible_Gene_Pairs"]:-NA}"
        LINE="$LINE,${row_map["Association_Rate"]:-NA}"
        LINE="$LINE,${row_map["Avg_Degree"]:-NA}"
        LINE="$LINE,${row_map["Module_Count"]:-NA}"
        LINE="$LINE,${row_map["Avg_Genes_per_Module"]:-NA}"
        LINE="$LINE,${row_map["Modularity"]:-NA}"
        ;;
      "Goldfinder")
        LINE="$LINE,${row_map["Dataset"]:-NA}"
        LINE="$LINE,${row_map["Total_Pangenome_Genes"]:-NA}"
        LINE="$LINE,${row_map["Network_Nodes"]:-NA}"
        LINE="$LINE,${row_map["Network_Edges"]:-NA}"
        LINE="$LINE,${row_map["Gene_Associations"]:-NA}"
        LINE="$LINE,NA" # Possible_Gene_Pairs
        LINE="$LINE,NA" # Association_Rate
        LINE="$LINE,NA" # Avg_Degree
        LINE="$LINE,${row_map["Cluster_Count"]:-NA}"       # Module_Count
        LINE="$LINE,${row_map["Avg_Cluster_Size"]:-NA}"    # Avg_Genes_per_Module
        LINE="$LINE,NA" # Modularity
        ;;
      *)
        echo "[⚠️] Unknown tool: $TOOL"
        continue
        ;;
    esac

    echo "$LINE" >> combined_summary.csv
  done
}

combine_by_header "Coinfinder" coinfinder_summary.csv
combine_by_header "Goldfinder" goldfinder_summary.csv

echo "[✅] Combined summary written to combined_summary.csv"
