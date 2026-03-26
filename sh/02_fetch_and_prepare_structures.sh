#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <pdb_token_list> <out_dir> [min_length] [dedup: yes|no]"
  exit 1
fi

TOKEN_LIST="$1"
OUT_DIR="$2"
MIN_LENGTH="${3:-0}"
DEDUP="${4:-yes}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p "$OUT_DIR"

RAW_FASTA="$OUT_DIR/structures.raw.fasta"
FAILED_TXT="$OUT_DIR/structures.failed.txt"
CLEAN_FASTA="$OUT_DIR/structures.cleaned.fasta"
CLEAN_TSV="$OUT_DIR/structures.cleaned.tsv"

python3 "$PROJECT_ROOT/scripts/core/getpdb.py" \
  "$TOKEN_LIST" \
  -o "$RAW_FASTA" \
  --failed "$FAILED_TXT"

DEDUP_FLAG="--dedup-by-sequence"
if [[ "$DEDUP" == "no" ]]; then
  DEDUP_FLAG="--no-dedup-by-sequence"
fi

python3 "$PROJECT_ROOT/scripts/core/04_prepare_structure_fasta.py" \
  "$RAW_FASTA" \
  -o "$CLEAN_FASTA" \
  --min-length "$MIN_LENGTH" \
  "$DEDUP_FLAG" \
  --report "$CLEAN_TSV"
