#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cog_short_fasta> <structure_fasta> <result_dir>"
  exit 1
fi

COG_FASTA="$1"
STRUCT_FASTA="$2"
RESULT_DIR="$3"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MUSCLE_BIN="${MUSCLE_BIN:-muscle}"
STRICT_TARGET_TOTAL="${STRICT_TARGET_TOTAL:-500}"
STRICT_THRESHOLD="${STRICT_THRESHOLD:-0.22}"
CORE_OCCUPANCY="${CORE_OCCUPANCY:-0.5}"
IQTREE_BIN="${IQTREE_BIN:-iqtree3}"
TREE_MODEL="${TREE_MODEL:-LG+I+R5}"

if ! command -v "$IQTREE_BIN" >/dev/null 2>&1; then
  if [[ -x /tmp/iqtree-env/bin/iqtree3 ]]; then
    IQTREE_BIN="/tmp/iqtree-env/bin/iqtree3"
  elif command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN="iqtree2"
  else
    echo "IQ-TREE binary was not found."
    echo "Set IQTREE_BIN explicitly or install iqtree3/iqtree2 in PATH."
    exit 1
  fi
fi

mkdir -p "$RESULT_DIR/aux"

python3 "$PROJECT_ROOT/scripts/core/04_prepare_structure_fasta.py" \
  "$STRUCT_FASTA" \
  -o "$RESULT_DIR/structures.cleaned.fasta" \
  --min-length 400 \
  --report "$RESULT_DIR/structures.cleaned.tsv"

cat "$COG_FASTA" "$RESULT_DIR/structures.cleaned.fasta" > "$RESULT_DIR/mega.input.fasta"

"$MUSCLE_BIN" -align "$RESULT_DIR/mega.input.fasta" -output "$RESULT_DIR/mega.msa.fasta"

python3 "$PROJECT_ROOT/scripts/core/06_filter_msa_keep_all_cogs.py" \
  "$RESULT_DIR/mega.msa.fasta" \
  --cog-fasta "$COG_FASTA" \
  -o "$RESULT_DIR/mega.filtered_keepcogs.fasta" \
  --threshold "$STRICT_THRESHOLD" \
  --target-total "$STRICT_TARGET_TOTAL" \
  --report "$RESULT_DIR/mega.filter.tsv"

python3 "$PROJECT_ROOT/scripts/core/14_filter_msa_columns_by_gap.py" \
  "$RESULT_DIR/mega.filtered_keepcogs.fasta" \
  -o "$RESULT_DIR/mega.core.msa.fasta" \
  --min-occupancy "$CORE_OCCUPANCY" \
  --report "$RESULT_DIR/mega.core.report.txt"

python3 "$PROJECT_ROOT/scripts/core/09_build_iqtree_tree.py" \
  "$RESULT_DIR/mega.core.msa.fasta" \
  -o "$RESULT_DIR/mega.tree.nwk" \
  --iqtree-bin "$IQTREE_BIN" \
  --model "$TREE_MODEL" \
  --bootstrap 0 \
  --alrt 0 \
  --threads 1 \
  --fast \
  --prefix "$RESULT_DIR/aux/mega"

python3 "$PROJECT_ROOT/scripts/core/02_sort_fasta_by_newick.py" \
  "$RESULT_DIR/mega.filtered_keepcogs.fasta" \
  "$RESULT_DIR/mega.tree.nwk" \
  -o "$RESULT_DIR/mega.reordered.fasta"

echo "Strict keep-all-COGs workflow done: $RESULT_DIR"
