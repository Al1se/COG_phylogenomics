#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cog_fasta> <cog_id> <result_dir>"
  exit 1
fi

COG_FASTA="$1"
COG_ID="$2"
RESULT_DIR="$3"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFERENCE_QUERY_COUNT="${REFERENCE_QUERY_COUNT:-3}"
REFERENCE_MOTIF="${REFERENCE_MOTIF:-}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
PDB_SEARCH_EVALUE="${PDB_SEARCH_EVALUE:-1e-10}"
PDB_SEARCH_IDENTITY="${PDB_SEARCH_IDENTITY:-0.30}"
PDB_SEARCH_MAX_HITS="${PDB_SEARCH_MAX_HITS:-100}"

mkdir -p "$RESULT_DIR"

QUERY_FASTA="$RESULT_DIR/reference_queries.fasta"
QUERY_REPORT="$RESULT_DIR/reference_queries.tsv"
TOKEN_LIST="$RESULT_DIR/pdb_tokens.similarity.txt"
SEARCH_REPORT="$RESULT_DIR/pdb_sequence_hits.tsv"

if [[ -n "$REFERENCE_FASTA" && -f "$REFERENCE_FASTA" ]]; then
  cp "$REFERENCE_FASTA" "$QUERY_FASTA"
  printf "rank\theader\tlength\tcontains_motif\tdistance_from_median\n" > "$QUERY_REPORT"
else
  python3 "$PROJECT_ROOT/scripts/structure/select_reference_queries.py" \
    "$COG_FASTA" \
    -o "$QUERY_FASTA" \
    --report "$QUERY_REPORT" \
    --max-queries "$REFERENCE_QUERY_COUNT" \
    --motif "$REFERENCE_MOTIF"
fi

python3 "$PROJECT_ROOT/scripts/structure/search_rcsb_by_sequence.py" \
  "$QUERY_FASTA" \
  -o "$TOKEN_LIST" \
  --report "$SEARCH_REPORT" \
  --evalue-cutoff "$PDB_SEARCH_EVALUE" \
  --identity-cutoff "$PDB_SEARCH_IDENTITY" \
  --max-hits-per-query "$PDB_SEARCH_MAX_HITS"

echo "$TOKEN_LIST"
