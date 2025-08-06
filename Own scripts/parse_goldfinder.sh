#!/bin/bash
# Usage: ./parse_goldfinder.sh <goldfinder_dir> <name>

DIR="$1"
NAME="$2"

presence_file="simulation/${NAME}/gene_presence_absence.csv"
GENE_COUNT="NA"
ASSOCIATIONS="NA"
MEAN_BEFORE="NA"
MEAN_AFTER="NA"
CLUSTER_COUNT="NA"
AVG_CLUSTER_SIZE="NA"
NODES="NA"
EDGES="NA"

[ -f "$presence_file" ] && GENE_COUNT=$(tail -n +2 "$presence_file" | wc -l)
[ -f "$DIR/simultaneous_association_significant_pairs.csv" ] && \
  ASSOCIATIONS=$(tail -n +2 "$DIR/simultaneous_association_significant_pairs.csv" | wc -l)

[ -f "$DIR/score_distribution_before_additional_simulation.txt" ] && \
  MEAN_BEFORE=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_before_additional_simulation.txt")

[ -f "$DIR/score_distribution_after_additional_simulation.txt" ] && \
  MEAN_AFTER=$(awk '{ total += $1; count++ } END { if (count > 0) print total/count; else print "NA" }' "$DIR/score_distribution_after_additional_simulation.txt")

[ -f "$DIR/association_clusters.txt" ] && \
  CLUSTER_COUNT=$(grep -c '^>' "$DIR/association_clusters.txt") && \
  AVG_CLUSTER_SIZE=$(awk '/^>/ {if (size > 0) total += size; clusters++; size = 0} !/^>/ {size++} END { if (clusters > 0) print total/clusters; else print "NA" }' "$DIR/association_clusters.txt")

[ -f "$DIR/cytoscape_input.csv" ] && \
  EDGES=$(grep -v '^\s*$' "$DIR/cytoscape_input.csv" | tail -n +2 | wc -l) && \
  NODES=$(tail -n +2 "$DIR/cytoscape_input.csv" | cut -d',' -f1,2 | tr ',' '\n' | sort | uniq | wc -l)

echo "$GENE_COUNT,$ASSOCIATIONS,$MEAN_BEFORE,$MEAN_AFTER,$CLUSTER_COUNT,$AVG_CLUSTER_SIZE,$NODES,$EDGES"
