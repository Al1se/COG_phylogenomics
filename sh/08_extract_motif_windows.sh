#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <input_fasta> <motif> <output_fasta> [left_flank] [right_flank] [wildcard]"
  exit 1
fi

INPUT_FASTA="$1"
MOTIF="$2"
OUTPUT_FASTA="$3"
LEFT_FLANK="${4:-2}"
RIGHT_FLANK="${5:-2}"
WILDCARD="${6:-X}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python3 "$PROJECT_ROOT/scripts/core/03_extract_motif_windows.py" \
  "$INPUT_FASTA" \
  -o "$OUTPUT_FASTA" \
  --motif "$MOTIF" \
  --left "$LEFT_FLANK" \
  --right "$RIGHT_FLANK" \
  --wildcard "$WILDCARD"
