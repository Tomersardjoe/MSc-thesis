#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i GPA_ROOT -o SUMMARIES_DIR -c COINFINDER_ROOT -g GOLDFINDER_ROOT -p PANFOREST_ROOT [-f]

Required:
  -i   Path to GPA root (e.g., simulation)
  -o   Path to output dir (e.g., summaries)
  -c   Path to coinfinder root (e.g., coinfinder)
  -g   Path to goldfinder root (e.g., goldfinder/goldfinder)
  -p   Path to panforest root (e.g., panforest)

Optional:
  -f   Force re-run even if SUMMARIES_DIR/{dataset}_metrics.csv exists
EOF
  exit 1
}

FORCE=0
GPA_ROOT=""
SUMMARIES_DIR=""
COINFINDER_ROOT=""
GOLDFINDER_ROOT=""
PANFOREST_ROOT=""

while getopts "i:o:c:g:p:f" opt; do
  case "$opt" in
    i) GPA_ROOT=$OPTARG ;;
    o) SUMMARIES_DIR=$OPTARG ;;
    c) COINFINDER_ROOT=$OPTARG ;;
    g) GOLDFINDER_ROOT=$OPTARG ;;
    p) PANFOREST_ROOT=$OPTARG ;;
    f) FORCE=1 ;;
    *) usage ;;
  esac
done

# Validate required args
if [[ -z "${GPA_ROOT}" || -z "${SUMMARIES_DIR}" || -z "${COINFINDER_ROOT}" || -z "${GOLDFINDER_ROOT}" || -z "${PANFOREST_ROOT}" ]]; then
  usage
fi

# Ensure output dir exists
mkdir -p "${SUMMARIES_DIR}"

# Discover dataset names: immediate subdirectories of GPA_ROOT (ignores hidden)
shopt -s nullglob
datasets=()
for d in "${GPA_ROOT}"/*/; do
  [[ -d "$d" ]] || continue
  datasets+=("$(basename "${d%/}")")
done
shopt -u nullglob

# Sort for deterministic ordering
if ((${#datasets[@]} == 0)); then
  echo "No dataset directories found under '${GPA_ROOT}'."
  exit 0
fi
IFS=$'\n' read -r -d '' -a datasets < <(printf '%s\n' "${datasets[@]}" | LC_ALL=C sort && printf '\0')
unset IFS

printf 'Datasets found:\n'
printf '  %s\n' "${datasets[@]}"

echo "Outputs will be written to '${SUMMARIES_DIR}'."

processed=0
skipped=0
failed=()

for dataset in "${datasets[@]}"; do
  summary_csv="${SUMMARIES_DIR}/${dataset}_metrics.csv"

  if [[ -f "${summary_csv}" && "${FORCE}" -eq 0 ]]; then
    echo "[SKIP] ${dataset}"
    ((++skipped))
    continue
  fi

  echo ">>> Processing: ${dataset}"

  (
    Rscript scripts/network_analysis.R \
      "${GPA_ROOT}" "${SUMMARIES_DIR}" \
      "${COINFINDER_ROOT}" "${GOLDFINDER_ROOT}" \
      "${PANFOREST_ROOT}" "${dataset}"
  )
  status=$?

  if [[ $status -ne 0 ]]; then
    echo "[FAIL] ${dataset} (exit ${status})"
    failed+=("${dataset}")
    continue
  fi
  ((++processed))
done

echo "-----"
echo "Processed: ${processed}"
echo "Skipped:   ${skipped}"
if ((${#failed[@]} > 0)); then
  echo "Failed:    ${#failed[@]} -> ${failed[*]}"
  exit 1
fi
echo "All requested datasets completed."

