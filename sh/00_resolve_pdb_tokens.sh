#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <COG_ID> <output_token_list> [mapping_file] [ligands_table]"
  exit 1
fi

COG_ID="$1"
OUTPUT_LIST="$2"
MAPPING_FILE="${3:-}"
LIGANDS_TABLE="${4:-}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "$MAPPING_FILE" ]]; then
  MAPPING_FILE="$PROJECT_ROOT/resources/supervisor/COG_to_PDB_chains.txt"
fi

if [[ ! -f "$MAPPING_FILE" ]]; then
  echo "Mapping file not found: $MAPPING_FILE"
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_LIST")"

CMD=(
  python3
  "$PROJECT_ROOT/scripts/supervisor/extract_cog_pdb_tokens.py"
  "$COG_ID"
  -m "$MAPPING_FILE"
  -o "$OUTPUT_LIST"
)

if [[ -n "$LIGANDS_TABLE" && -f "$LIGANDS_TABLE" ]]; then
  CMD+=(--ligands-table "$LIGANDS_TABLE")
fi

"${CMD[@]}"
