#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <token_list> <mmcif_dir> <contacts_tsv> <reordered_msa> <out_dir> [sites_tsv] [query_fasta] [title]"
  exit 1
fi

TOKEN_LIST="$1"
MMCIF_DIR="$2"
CONTACTS_TSV="$3"
REORDERED_MSA="$4"
OUT_DIR="$5"
SITES_TSV="${6:-}"
QUERY_FASTA="${7:-}"
TITLE="${8:-STRUCTURE SITE COORDINATES IN REORDERED MSA}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p "$OUT_DIR"

TOKEN_MAP="$OUT_DIR/structure_sequence_mapping.tsv"
MAPPED_CONTACTS="$OUT_DIR/contacts_in_reordered.tsv"
SITE_COORDS_TXT="$OUT_DIR/site_coordinates_in_reordered.txt"
SITE_COORDS_TSV="$OUT_DIR/site_coordinates_in_reordered.tsv"

python3 "$PROJECT_ROOT/scripts/structure/map_structure_sequences_to_alignment.py" \
  "$TOKEN_LIST" \
  --mmcif-dir "$MMCIF_DIR" \
  --msa-fasta "$REORDERED_MSA" \
  -o "$TOKEN_MAP"

python3 "$PROJECT_ROOT/scripts/structure/map_contacts_via_sequence_matches.py" \
  "$CONTACTS_TSV" \
  "$TOKEN_MAP" \
  "$REORDERED_MSA" \
  -o "$MAPPED_CONTACTS"

SUMMARY_ARGS=(
  "$PROJECT_ROOT/scripts/structure/summarize_structure_sites_in_alignment.py"
  "$MAPPED_CONTACTS"
  --msa-fasta "$REORDERED_MSA"
  --output-txt "$SITE_COORDS_TXT"
  --output-tsv "$SITE_COORDS_TSV"
  --title "$TITLE"
)

if [[ -n "$SITES_TSV" && -f "$SITES_TSV" ]]; then
  SUMMARY_ARGS+=(--sites-tsv "$SITES_TSV")
fi

if [[ -n "$QUERY_FASTA" && -f "$QUERY_FASTA" ]]; then
  SUMMARY_ARGS+=(--query-fasta "$QUERY_FASTA")
fi

python3 "${SUMMARY_ARGS[@]}"
