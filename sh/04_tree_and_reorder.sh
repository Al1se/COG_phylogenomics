#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <combined_msa> <out_dir> [core_occupancy]"
  exit 1
fi

COMBINED_MSA="$1"
OUT_DIR="$2"
CORE_OCCUPANCY="${3:-0.7}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
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

mkdir -p "$OUT_DIR/aux"

CORE_MSA="$OUT_DIR/combined.core.msa.fasta"
CORE_REPORT="$OUT_DIR/combined.core.report.txt"
TREE_NWK="$OUT_DIR/tree.nwk"
REORDERED="$OUT_DIR/reordered.fasta"

python3 "$PROJECT_ROOT/scripts/core/14_filter_msa_columns_by_gap.py" \
  "$COMBINED_MSA" \
  -o "$CORE_MSA" \
  --min-occupancy "$CORE_OCCUPANCY" \
  --report "$CORE_REPORT"

python3 "$PROJECT_ROOT/scripts/core/09_build_iqtree_tree.py" \
  "$CORE_MSA" \
  -o "$TREE_NWK" \
  --iqtree-bin "$IQTREE_BIN" \
  --model "$TREE_MODEL" \
  --bootstrap 0 \
  --alrt 0 \
  --threads 1 \
  --fast \
  --prefix "$OUT_DIR/aux/finaltree"

python3 "$PROJECT_ROOT/scripts/core/02_sort_fasta_by_newick.py" \
  "$COMBINED_MSA" \
  "$TREE_NWK" \
  -o "$REORDERED"
