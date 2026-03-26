#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cog_fasta> <structures_fasta> <out_dir>"
  exit 1
fi

COG_FASTA="$1"
STRUCT_FASTA="$2"
OUT_DIR="$3"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MUSCLE_BIN="${MUSCLE_BIN:-muscle}"

mkdir -p "$OUT_DIR"

COMBINED_FASTA="$OUT_DIR/combined.input.fasta"
STRUCT_MSA="$OUT_DIR/structures.only.msa.fasta"
COMBINED_MSA="$OUT_DIR/combined.msa.fasta"

count_fasta_records() {
  local fasta_path="$1"
  if [[ ! -s "$fasta_path" ]]; then
    echo 0
    return
  fi
  grep -c '^>' "$fasta_path" || true
}

if [[ -s "$STRUCT_FASTA" ]]; then
  cat "$COG_FASTA" "$STRUCT_FASTA" > "$COMBINED_FASTA"
  STRUCT_COUNT="$(count_fasta_records "$STRUCT_FASTA")"
  if [[ "$STRUCT_COUNT" -ge 2 ]]; then
    "$MUSCLE_BIN" -align "$STRUCT_FASTA" -output "$STRUCT_MSA"
  else
    cp "$STRUCT_FASTA" "$STRUCT_MSA"
  fi
else
  cp "$COG_FASTA" "$COMBINED_FASTA"
  : > "$STRUCT_MSA"
fi

COMBINED_COUNT="$(count_fasta_records "$COMBINED_FASTA")"
if [[ "$COMBINED_COUNT" -ge 2 ]]; then
  "$MUSCLE_BIN" -align "$COMBINED_FASTA" -output "$COMBINED_MSA"
else
  cp "$COMBINED_FASTA" "$COMBINED_MSA"
fi
