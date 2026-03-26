#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <aligned_msa> <output_fasta> <threshold> <reference_id> [reference_id ...]"
  exit 1
fi

INPUT_MSA="$1"
OUTPUT_FASTA="$2"
THRESHOLD="$3"
shift 3

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ARGS=()
for ref in "$@"; do
  ARGS+=(--reference-id "$ref")
done

python3 "$PROJECT_ROOT/scripts/core/01_filter_msa_by_similarity.py" \
  "$INPUT_MSA" \
  -o "$OUTPUT_FASTA" \
  --threshold "$THRESHOLD" \
  "${ARGS[@]}"
