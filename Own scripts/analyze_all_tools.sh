#!/bin/bash

echo "[🚀] Starting tool analyses..."

### 📊 Parse Coinfinder Output
echo "[📦] Parsing Coinfinder results..."
echo "Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity" > coinfinder_summary.csv

for DIR in ./coinfinder/*/; do
  NAME=$(basename "$DIR")
  if METRICS=$(./parse_coinfinder.sh "$DIR" "$NAME"); then
    echo "$NAME,$METRICS" >> coinfinder_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> coinfinder_summary.csv
    echo "[⚠] Failed to parse Coinfinder: $NAME — see debug output for details"
  fi
  echo "[✔] Finished parsing Coinfinder: $NAME"
done

### 🧬 Parse Goldfinder Output (Associations only)
echo "[📦] Parsing Goldfinder results..."
echo "Dataset,Gene_Associations,Mean_Score_Before,Mean_Score_After,Cluster_Count,Avg_Cluster_Size,Network_Nodes,Network_Edges" > goldfinder_summary.csv

for DIR in ./goldfinder/goldfinder/*/; do
  NAME=$(basename "$DIR")
  ASSOCIATIONS="NA"
  MEAN_BEFORE="NA"
  MEAN_AFTER="NA"
  CLUSTER_COUNT="NA"
  AVG_CLUSTER_SIZE="NA"
  NODES="NA"
  EDGES="NA"

  [ -f "$DIR/simultaneous_association_significant_pairs.csv" ] && \
    ASSOCIATIONS=$(tail -n +2 "$DIR/simultaneous_association_significant_pairs.csv" | wc -l)

  [ -f "$DIR/score_distribution_before_additional_simulation.txt" ] && \
    MEAN_BEFORE=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_before_additional_simulation.txt")

  [ -f "$DIR/score_distribution_after_additional_simulation.txt" ] && \
    MEAN_AFTER=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_after_additional_simulation.txt")

  [ -f "$DIR/association_clusters.txt" ] && \
    CLUSTER_COUNT=$(grep -c '^>' "$DIR/association_clusters.txt") && \
    AVG_CLUSTER_SIZE=$(awk '/^>/ {if (size > 0) total += size; clusters++ ; size = 0 } !/^>/ {size++} END { if (clusters > 0) print total/clusters; else print "NA" }' "$DIR/association_clusters.txt")

  [ -f "$DIR/cytoscape_input.csv" ] && \
    EDGES=$(tail -n +2 "$DIR/cytoscape_input.csv" | wc -l) && \
    NODES=$(tail -n +2 "$DIR/cytoscape_input.csv" | cut -d',' -f1,2 | tr ',' '\n' | sort | uniq | wc -l)

  echo "$NAME,$ASSOCIATIONS,$MEAN_BEFORE,$MEAN_AFTER,$CLUSTER_COUNT,$AVG_CLUSTER_SIZE,$NODES,$EDGES" >> goldfinder_summary.csv
  echo "[✔] Finished parsing Goldfinder: $NAME"
done

### 🧺 Combine Results into a Unified Summary
echo "[🧩] Combining tool summaries..."

echo "Tool,Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity,Gene_Associations,Mean_Score_Before,Mean_Score_After,Cluster_Count,Avg_Cluster_Size,Network_Nodes,Network_Edges" > combined_summary.csv

tail -n +2 coinfinder_summary.csv | while IFS= read -r line; do
  echo "Coinfinder,$line,NA,NA,NA,NA,NA,NA,NA" >> combined_summary.csv
done

tail -n +2 goldfinder_summary.csv | while IFS= read -r line; do
  echo "$line" | awk -F',' '{printf "Goldfinder,%s", $0; for(i=NF+1;i<=17;i++) printf ",NA"; print ""}' >> combined_summary.csv
done

echo "[✅] Combined summary written to combined_summary.csv"
