#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR"/.. && pwd)"
COINFINDER_DIR="$PROJECT_ROOT/coinfinder"
SIMULATION_DIR="$PROJECT_ROOT/simulation"
DVAL_DIR="$PROJECT_ROOT/d_values"

echo "[🚀] Checking for missing d_values.tsv files..."

mkdir -p "$DVAL_DIR"

for DIR in "$COINFINDER_DIR"/*/; do
  NAME=$(basename "$DIR")
  DVAL_FILE="$DVAL_DIR/${NAME}_d_values.tsv"

  # Skip underscore-prefixed dirs like __pycache__
  [[ "$NAME" =~ ^_+ ]] && echo "[⏭] Skipping $NAME" && continue

  if [[ -f "$DVAL_FILE" ]]; then
    echo "[✅] D-values already exist for: $NAME"
    continue
  fi

  ALPHA_BETA_FILE="$SIMULATION_DIR/$NAME/alpha_beta.tab"
  TREE_FILE="$SIMULATION_DIR/$NAME/tree.nwk"

  if [[ ! -f "$ALPHA_BETA_FILE" ]] || [[ ! -f "$TREE_FILE" ]]; then
    echo "[⚠] Missing input files for: $NAME"
    continue
  fi

  echo "[⚙️] Generating D-values for: $NAME"

  OUT_PREFIX="$DVAL_DIR/tmp_$NAME"
  DVAL_FILE="$DVAL_DIR/${NAME}_d_values.tsv"

  # Run Coinfinder
  coinfinder -i "$PRES_FILE" \
             -p "$TREE_FILE" \
             -x 4 -a -E \
             -o "$OUT_PREFIX"

  # Rename output file
  if [[ -f "${OUT_PREFIX}_nodes.tsv" ]]; then
    mv "${OUT_PREFIX}_nodes.tsv" "$DVAL_FILE"
    echo "[📄] Saved: $DVAL_FILE"
  else
    echo "[⚠] Coinfinder did not produce nodes.tsv for: $NAME"
  fi

  # Clean up other Coinfinder outputs
  rm -f "${OUT_PREFIX}_summary.tsv" "${OUT_PREFIX}_pairs.tsv" "${OUT_PREFIX}_edges.tsv" "${OUT_PREFIX}_gamma.tsv"
 
  echo "[✔] Finished D-value generation for: $NAME"
done

echo "[🎉] All missing d_values.tsv files generated!"
