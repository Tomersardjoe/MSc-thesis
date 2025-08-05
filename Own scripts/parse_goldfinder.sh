#!/bin/bash

# Script to parse Goldfinder association output and summarize metrics
# Output: goldfinder_summary.csv

echo "[🔍] Starting Goldfinder association parsing..."

# Output CSV header
echo "Dataset,Gene_Associations,Mean_Score_Before,Mean_Score_After,Cluster_Count,Avg_Cluster_Size,Network_Nodes,Network_Edges" > goldfinder_summary.csv

# Loop through datasets
for DIR in ./goldfinder/goldfinder/*; do
  if [ -d "$DIR" ]; then
    NAME=$(basename "$DIR")

    # Initialize metric variables
    ASSOCIATIONS="NA"
    MEAN_BEFORE="NA"
    MEAN_AFTER="NA"
    CLUSTER_COUNT="NA"
    AVG_CLUSTER_SIZE="NA"
    NODES="NA"
    EDGES="NA"

    # 1. Count associated gene pairs
    if [ -f "$DIR/simultaneous_association_significant_pairs.csv" ]; then
      ASSOCIATIONS=$(tail -n +2 "$DIR/simultaneous_association_significant_pairs.csv" | wc -l)
    fi

    # 2. Parse score distributions
    if [ -f "$DIR/score_distribution_before_additional_simulation.txt" ]; then
      MEAN_BEFORE=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_before_additional_simulation.txt")
    fi

    if [ -f "$DIR/score_distribution_after_additional_simulation.txt" ]; then
      MEAN_AFTER=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_after_additional_simulation.txt")
    fi

    # 3. Parse cluster count and average size
    if [ -f "$DIR/association_clusters.txt" ]; then
      CLUSTER_COUNT=$(grep -c '^>' "$DIR/association_clusters.txt")
      AVG_CLUSTER_SIZE=$(awk '/^>/ {if (size > 0) total += size; clusters++ ; size = 0 } !/^>/ {size++} END { if (clusters > 0) print total/clusters; else print "NA" }' "$DIR/association_clusters.txt")
    fi

    # 4. Cytoscape file: nodes and edges
    if [ -f "$DIR/cytoscape_input.csv" ]; then
      EDGES=$(tail -n +2 "$DIR/cytoscape_input.csv" | wc -l)
      NODES=$(tail -n +2 "$DIR/cytoscape_input.csv" | cut -d',' -f1,2 | tr ',' '\n' | sort | uniq | wc -l)
    fi

    # Save metrics
    echo "$NAME,$ASSOCIATIONS,$MEAN_BEFORE,$MEAN_AFTER,$CLUSTER_COUNT,$AVG_CLUSTER_SIZE,$NODES,$EDGES" >> goldfinder_summary.csv
    echo "[✔] Finished parsing $NAME"
  fi
done

echo "[🎯] All Goldfinder association results parsed. Summary saved to goldfinder_summary.csv"
