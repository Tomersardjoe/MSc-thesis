#!/bin/bash
echo "[🚀] Starting tool analyses..."

# Coinfinder
echo "[📦] Parsing Coinfinder results..."
echo "Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity" > coinfinder_summary.csv

for DIR in ./coinfinder/*/; do
  NAME=$(basename "$DIR")
  if METRICS=$(./parse_coinfinder.sh "$DIR" "$NAME"); then
    echo "$NAME,$METRICS" >> coinfinder_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> coinfinder_summary.csv
    echo "[⚠] Failed to parse Coinfinder: $NAME"
  fi
  echo "[✔] Finished parsing Coinfinder: $NAME"
done

# Goldfinder
echo "[📦] Parsing Goldfinder results..."
echo "Dataset,Total_Pangenome_Genes,Gene_Associations,Mean_Score_Before,Mean_Score_After,Cluster_Count,Avg_Cluster_Size,Network_Nodes,Network_Edges" > goldfinder_summary.csv

for DIR in ./goldfinder/goldfinder/*/; do
  NAME=$(basename "$DIR")
  [[ "$NAME" =~ ^_+ ]] && echo "[⏭] Skipping $NAME" && continue
  if METRICS=$(./parse_goldfinder.sh "$DIR" "$NAME"); then
    echo "$NAME,$METRICS" >> goldfinder_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA" >> goldfinder_summary.csv
    echo "[⚠] Failed to parse Goldfinder: $NAME"
  fi
  echo "[✔] Finished parsing Goldfinder: $NAME"
done

### 🔗 Combine Results Across Tools
echo "[🔁] Invoking combiner to generate combined summary..."
./combine_summaries.sh
