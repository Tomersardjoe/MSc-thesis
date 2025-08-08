#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Resolve the directory this script lives in and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR"/.. && pwd)"

# Usage check
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <simulation_dir> <dataset_prefix>" >&2
  exit 1
fi

SIM_DIR="$1"
PREFIX="$2"
DATASET="$PREFIX"

# File paths
PRESENCE_FILE="$PROJECT_ROOT/simulation/${PREFIX}/gene_presence_absence.csv"
PAIRS="$SIM_DIR/simultaneous_association_significant_pairs.csv"
CLUSTERS="$SIM_DIR/association_clusters.txt"
NETWORK="$SIM_DIR/cytoscape_input.csv"

# Validate inputs
if [[ ! -f "$PRESENCE_FILE" ]] \
  || [[ ! -f "$PAIRS" ]] \
  || [[ ! -f "$CLUSTERS" ]] \
  || [[ ! -f "$NETWORK" ]]; then
  echo "ERROR: Missing required files in $SIM_DIR or simulation/${PREFIX}" >&2
  exit 1
fi

# -------------------------------------------------------------------
# Total pangenome genes
# -------------------------------------------------------------------
TOTAL_PANGENOME_GENES=$(tail -n +2 "$PRESENCE_FILE" | wc -l)

# -------------------------------------------------------------------
# Network statistics
# -------------------------------------------------------------------
NETWORK_EDGES=$(grep -v '^\s*$' "$NETWORK" | tail -n +2 | wc -l)
NETWORK_NODES=$(tail -n +2 "$NETWORK" \
  | cut -d',' -f1,2 \
  | tr ',' '\n' \
  | sort -u \
  | wc -l)

# -------------------------------------------------------------------
# Association statistics
# -------------------------------------------------------------------
GENE_ASSOCIATIONS=$(tail -n +2 "$PAIRS" | wc -l)
POSSIBLE_PAIRS=$(( NETWORK_NODES * (NETWORK_NODES - 1) / 2 ))
ASSOCIATION_RATE=$(awk -v a="$GENE_ASSOCIATIONS" -v p="$POSSIBLE_PAIRS" \
  'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

# -------------------------------------------------------------------
# Extract unique associated genes with p-adj < 0.05
# -------------------------------------------------------------------
GOLDGENES="$SIM_DIR/goldfinder_genes.txt"
awk -F',' 'NR>1 && $4 < 0.05 {
  gsub(/"/, "", $1); gsub(/"/, "", $2);
  print $1; print $2;
}' "$PAIRS" | sort -u > "$GOLDGENES"

# -------------------------------------------------------------------
# Cluster statistics
# -------------------------------------------------------------------
MODULE_COUNT=$(grep -c '^>' "$CLUSTERS")

# -------------------------------------------------------------------
# Network metrics via Python
# -------------------------------------------------------------------
METRICS_OUTPUT=$(python3 "$SCRIPT_DIR"/network_metrics.py "$SIM_DIR")
AVG_DEGREE=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Avg\. Degree: \K[\d.]+')
MODULARITY=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Modularity: \K[\d.]+')
AVG_GENES_PER_MODULE=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Avg\. Module Size.*: \K[\d.]+')

# -------------------------------------------------------------------
# Emit a single CSV line (no header)
# -------------------------------------------------------------------
echo "${DATASET},${TOTAL_PANGENOME_GENES},${NETWORK_NODES},${POSSIBLE_PAIRS},${GENE_ASSOCIATIONS},${ASSOCIATION_RATE},${MODULE_COUNT},${AVG_GENES_PER_MODULE},${AVG_DEGREE},${MODULARITY}"
