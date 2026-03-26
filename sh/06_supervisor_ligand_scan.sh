#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <pdb_or_cif_dir> <output_prefix> [ligands_file] [chains_file]"
  exit 1
fi

PDB_DIR="$1"
OUTPUT_PREFIX="$2"
LIGANDS_FILE="${3:-}"
CHAINS_FILE="${4:-}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUPERVISOR_PYTHON_BIN="${SUPERVISOR_PYTHON_BIN:-python2}"

if [[ -z "$LIGANDS_FILE" ]]; then
  LIGANDS_FILE="$PROJECT_ROOT/resources/supervisor/req_ligands.txt"
fi

if [[ -z "$CHAINS_FILE" ]]; then
  CHAINS_FILE="$PROJECT_ROOT/resources/supervisor/COG_to_PDB_chains.txt"
fi

if ! command -v "$SUPERVISOR_PYTHON_BIN" >/dev/null 2>&1; then
  echo "Supervisor tools require a Python 2 interpreter."
  echo "Set SUPERVISOR_PYTHON_BIN to a valid executable, for example python2.7."
  exit 1
fi

"$SUPERVISOR_PYTHON_BIN" "$PROJECT_ROOT/scripts/supervisor/find_ligands_in_PDB.py" \
  -i "$PDB_DIR" \
  -o "$OUTPUT_PREFIX" \
  -l "$LIGANDS_FILE" \
  -c "$CHAINS_FILE"
