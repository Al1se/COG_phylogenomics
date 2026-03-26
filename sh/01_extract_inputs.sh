#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <COG_ID> <out_dir> <mode: full|short|representative>"
  exit 1
fi

COG_ID="$1"
OUT_DIR="$2"
MODE="$3"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/input}"
SHORT_LIST="${SHORT_LIST:-$PROJECT_ROOT/resources/shortlists/275list.txt}"
MOTIF="${MOTIF:-NADFDGD}"

mkdir -p "$OUT_DIR"

case "$MODE" in
  full)
    python3 "$PROJECT_ROOT/scripts/core/cog-extractor-full.py" \
      "$COG_ID" \
      --cog-csv "$DATA_DIR/cog-24.cog.csv" \
      --proteins-fasta "$DATA_DIR/2296Genomes.prot.fasta" \
      -o "$OUT_DIR/${COG_ID}.full.raw.fasta"
    ;;
  short)
    python3 "$PROJECT_ROOT/scripts/core/cog-extractor-short-list.py" \
      "$COG_ID" \
      --short-list "$SHORT_LIST" \
      --cog-csv "$DATA_DIR/cog-24.cog.csv" \
      --proteins-fasta "$DATA_DIR/2296Genomes.prot.fasta" \
      -o "$OUT_DIR/${COG_ID}.short.raw.fasta"
    ;;
  representative)
    python3 "$PROJECT_ROOT/scripts/core/12_extract_cog_representatives.py" \
      "$COG_ID" \
      --short-list "$SHORT_LIST" \
      --cog-csv "$DATA_DIR/cog-24.cog.csv" \
      --proteins-fasta "$DATA_DIR/2296Genomes.prot.fasta" \
      --motif "$MOTIF" \
      -o "$OUT_DIR/${COG_ID}.representatives.fasta" \
      --report "$OUT_DIR/${COG_ID}.representatives.tsv"
    ;;
  *)
    echo "Unknown mode: $MODE"
    exit 1
    ;;
esac
