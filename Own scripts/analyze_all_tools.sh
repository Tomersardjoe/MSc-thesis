#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Resolve directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR"/.. && pwd)"
SUMMARIES_DIR="$PROJECT_ROOT/summaries"

# Create summaries directory if it doesn't exist already
mkdir -p "$SUMMARIES_DIR"

echo "[🚀] Starting tool analyses..."

# -----------------------------------------------------------------------------
# Coinfinder
# -----------------------------------------------------------------------------
echo "[📦] Parsing Coinfinder results..."
echo "Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity" \
  > "$SUMMARIES_DIR"/coinfinder_summary.csv

for DIR in "$PROJECT_ROOT"/coinfinder/*/; do
  NAME=$(basename "$DIR")

  if METRICS=$("$SCRIPT_DIR"/parse_coinfinder.sh "$DIR" "$NAME"); then
    echo "$METRICS" >> "$SUMMARIES_DIR"/coinfinder_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" \
      >> "$SUMMARIES_DIR"/coinfinder_summary.csv
    echo "[⚠] Failed to parse Coinfinder: $NAME"
  fi

  echo "[✔] Finished parsing Coinfinder: $NAME"
done

# -----------------------------------------------------------------------------
# Goldfinder
# -----------------------------------------------------------------------------
echo "[📦] Parsing Goldfinder results..."
echo "Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity" \
  > "$SUMMARIES_DIR"/goldfinder_summary.csv

for DIR in "$PROJECT_ROOT"/goldfinder/goldfinder/*/; do
  NAME=$(basename "$DIR")
  [[ "$NAME" =~ ^_+ ]] && echo "[⏭] Skipping $NAME" && continue

  if METRICS=$("$SCRIPT_DIR"/parse_goldfinder.sh "$DIR" "$NAME"); then
    echo "$METRICS" >> "$SUMMARIES_DIR"/goldfinder_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" \
      >> "$SUMMARIES_DIR"/goldfinder_summary.csv
    echo "[⚠] Failed to parse Goldfinder: $NAME"
  fi

  echo "[✔] Finished parsing Goldfinder: $NAME"
done

# -----------------------------------------------------------------------------
# Panforest
# -----------------------------------------------------------------------------
echo "[📦] Parsing Panforest results..."
echo "Dataset,Total_Pangenome_Genes,Total_Networked_Genes,Possible_Gene_Pairs,Significant_Associations,Association_Rate,Module_Count,Avg_Genes_per_Module,Avg_Degree,Modularity" \
  > "$SUMMARIES_DIR"/panforest_summary.csv

for DIR in "$PROJECT_ROOT"/panforest/*/; do
  NAME=$(basename "$DIR")
  [[ "$NAME" =~ ^_+ ]] && echo "[⏭] Skipping $NAME" && continue
  NAME=$(basename "$DIR")
  PRES_FILE="$PROJECT_ROOT/simulation/$NAME/gene_presence_absence.csv"
  IMP_FILE="$DIR/imp.csv"
  GRAPHML_OUT="$DIR/filtered_network.graphml"

  if [[ ! -f "$PRES_FILE" ]] || [[ ! -f "$IMP_FILE" ]]; then
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$SUMMARIES_DIR"/panforest_summary.csv
    echo "[⚠] Missing required files for Panforest: $NAME"
    continue
  fi

  # Capture metrics line printed by R script and append to summary
  if METRICS=$(Rscript "$SCRIPT_DIR"/parse_panforest.R "$PRES_FILE" "$IMP_FILE" "$GRAPHML_OUT"); then
    echo "$METRICS" >> "$SUMMARIES_DIR"/panforest_summary.csv
  else
    echo "$NAME,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$SUMMARIES_DIR"/panforest_summary.csv
    echo "[⚠] Failed Panforest analysis: $NAME"
  fi

  echo "[✔] Finished parsing Panforest: $NAME"
done

# -----------------------------------------------------------------------------
# Combine Summaries
# -----------------------------------------------------------------------------
echo "[🔁] Invoking combiner to generate combined summary..."
bash "$SCRIPT_DIR"/combine_summaries.sh \
  "$SUMMARIES_DIR"/coinfinder_summary.csv \
  "$SUMMARIES_DIR"/goldfinder_summary.csv \
  "$SUMMARIES_DIR"/panforest_summary.csv \
  "$SUMMARIES_DIR"/combined_summary.csv

echo "[🎉] All analyses complete!"
