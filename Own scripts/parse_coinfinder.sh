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
NODES="${SIM_DIR}/${PREFIX}_nodes.tsv"
PAIRS="${SIM_DIR}/${PREFIX}_pairs.tsv"
COMPONENTS="${SIM_DIR}/${PREFIX}_components.tsv"
PRESENCE_FILE="${PROJECT_ROOT}/simulation/${PREFIX}/gene_presence_absence.csv"

# Validate inputs
if [[ ! -f "$NODES" ]] || [[ ! -f "$PAIRS" ]] || [[ ! -f "$COMPONENTS" ]] || [[ ! -f "$PRESENCE_FILE" ]]; then
  echo "ERROR: Missing required files in $SIM_DIR or simulation/${PREFIX}" >&2
  exit 1
fi

# Total pangenome genes
TOTAL_PANGENOME_GENES=$(tail -n +2 "$PRESENCE_FILE" | wc -l)

# Networked genes
NETWORK_NODES=$(tail -n +2 "$NODES" \
  | cut -f1 \
  | sort -u \
  | wc -l)

# Association statistics
GENE_ASSOCIATIONS=$(awk '$3 < 0.05' "$PAIRS" | wc -l)
POSSIBLE_PAIRS=$(( NETWORK_NODES * (NETWORK_NODES - 1) / 2 ))
ASSOCIATION_RATE=$(awk -v a="$GENE_ASSOCIATIONS" -v p="$POSSIBLE_PAIRS" \
  'BEGIN { if (p > 0) printf("%.2f", (a/p)*100); else print "0.00" }')

# Module statistics
MODULE_COUNT=$(cut -f2 "$COMPONENTS" | sort -u | wc -l)

# Gather additional metrics via Python
METRICS_OUTPUT=$(python3 "$SCRIPT_DIR"/network_metrics.py "$SIM_DIR")
AVG_GENES_PER_MODULE=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Avg\. Module Size.*: \K[\d.]+')
AVG_DEGREE=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Avg\. Degree: \K[\d.]+')
MODULARITY=$(echo "$METRICS_OUTPUT" \
  | grep -oP 'Modularity: \K[\d.]+')

# Emit a single CSV line (no header)
echo "${DATASET},${TOTAL_PANGENOME_GENES},${NETWORK_NODES},${POSSIBLE_PAIRS},${GENE_ASSOCIATIONS},${ASSOCIATION_RATE},${MODULE_COUNT},${AVG_GENES_PER_MODULE},${AVG_DEGREE},${MODULARITY}"
