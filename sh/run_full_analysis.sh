#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <COG_ID> <result_dir> <mode: full|short|representative> [pdb_token_list]"
  exit 1
fi

COG_ID="$1"
RESULT_DIR="$2"
MODE="$3"
PDB_LIST="${4:-}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STRUCT_MIN_LENGTH="${STRUCT_MIN_LENGTH:-0}"
STRUCT_DEDUP="${STRUCT_DEDUP:-yes}"
CORE_OCCUPANCY="${CORE_OCCUPANCY:-0.7}"
RUN_STRUCTURE_MAPPING="${RUN_STRUCTURE_MAPPING:-1}"
RUN_SUPERVISOR_SCAN="${RUN_SUPERVISOR_SCAN:-auto}"
METALS="${METALS:-ZN,MG}"
SUPERVISOR_PYTHON_BIN="${SUPERVISOR_PYTHON_BIN:-python2}"

mkdir -p "$RESULT_DIR"

bash "$PROJECT_ROOT/sh/01_extract_inputs.sh" "$COG_ID" "$RESULT_DIR" "$MODE"

case "$MODE" in
  full)
    COG_FASTA="$RESULT_DIR/${COG_ID}.full.raw.fasta"
    ;;
  short)
    COG_FASTA="$RESULT_DIR/${COG_ID}.short.raw.fasta"
    ;;
  representative)
    COG_FASTA="$RESULT_DIR/${COG_ID}.representatives.fasta"
    ;;
  *)
    echo "Unknown mode: $MODE"
    exit 1
    ;;
esac

COG_SUFFIX="${COG_ID#COG}"
PDB_SOURCE="none"

if [[ -n "$PDB_LIST" && -f "$PDB_LIST" ]]; then
  PDB_SOURCE="explicit"
fi

if [[ -z "$PDB_LIST" ]]; then
  DISCOVERED_PDB_LIST="$RESULT_DIR/pdb_tokens.similarity.txt"
  if bash "$PROJECT_ROOT/sh/00_discover_pdb_tokens.sh" "$COG_FASTA" "$COG_ID" "$RESULT_DIR" >/dev/null 2>&1; then
    if [[ -s "$DISCOVERED_PDB_LIST" ]]; then
      PDB_LIST="$DISCOVERED_PDB_LIST"
      PDB_SOURCE="sequence_similarity"
    fi
  fi
fi

if [[ -z "$PDB_LIST" ]]; then
  for candidate in \
    "$PROJECT_ROOT/resources/pdb_lists/$COG_ID" \
    "$PROJECT_ROOT/resources/pdb_lists/ZN$COG_SUFFIX" \
    "$PROJECT_ROOT/resources/pdb_lists/MG$COG_SUFFIX"
  do
    if [[ -f "$candidate" ]]; then
      PDB_LIST="$candidate"
      PDB_SOURCE="resources"
      break
    fi
  done
fi

if [[ -z "$PDB_LIST" ]]; then
  AUTO_PDB_LIST="$RESULT_DIR/pdb_tokens.auto.txt"
  bash "$PROJECT_ROOT/sh/00_resolve_pdb_tokens.sh" "$COG_ID" "$AUTO_PDB_LIST"
  if [[ -s "$AUTO_PDB_LIST" ]]; then
    PDB_LIST="$AUTO_PDB_LIST"
    PDB_SOURCE="supervisor_mapping"
  fi
fi

printf "%s\n" "$PDB_SOURCE" > "$RESULT_DIR/pdb_token_source.txt"

if [[ -n "$PDB_LIST" && -f "$PDB_LIST" ]]; then
  echo "Using PDB token list: $PDB_LIST"
  echo "PDB token source: $PDB_SOURCE"
  bash "$PROJECT_ROOT/sh/02_fetch_and_prepare_structures.sh" "$PDB_LIST" "$RESULT_DIR" "$STRUCT_MIN_LENGTH" "$STRUCT_DEDUP"
else
  echo "No PDB token list was provided or auto-detected; running without structure sequences."
  : > "$RESULT_DIR/structures.raw.fasta"
  : > "$RESULT_DIR/structures.cleaned.fasta"
  : > "$RESULT_DIR/structures.failed.txt"
fi

bash "$PROJECT_ROOT/sh/03_build_msa.sh" "$COG_FASTA" "$RESULT_DIR/structures.cleaned.fasta" "$RESULT_DIR"
bash "$PROJECT_ROOT/sh/04_tree_and_reorder.sh" "$RESULT_DIR/combined.msa.fasta" "$RESULT_DIR" "$CORE_OCCUPANCY"

if [[ "$RUN_STRUCTURE_MAPPING" == "1" && -n "$PDB_LIST" && -f "$PDB_LIST" ]]; then
  MMCIF_DIR="$RESULT_DIR/mmcif"
  PDB_DIR="$RESULT_DIR/pdb"
  PYMOL_DIR="$RESULT_DIR/pymol"
  SITE_TITLE="${COG_ID} STRUCTURE SITE COORDINATES IN REORDERED MSA"
  python3 "$PROJECT_ROOT/scripts/structure/download_rcsb_mmcif.py" "$PDB_LIST" -o "$MMCIF_DIR"
  python3 "$PROJECT_ROOT/scripts/structure/download_rcsb_pdb.py" "$PDB_LIST" -o "$PDB_DIR" || true
  python3 "$PROJECT_ROOT/scripts/structure/analyze_metal_contacts.py" \
    "$PDB_LIST" \
    --mmcif-dir "$MMCIF_DIR" \
    --metals "$METALS" \
    --contacts-tsv "$RESULT_DIR/structure_contacts.tsv" \
    --sites-tsv "$RESULT_DIR/structure_sites.tsv" \
    --pymol-dir "$PYMOL_DIR"
  bash "$PROJECT_ROOT/sh/05_map_structure_sites.sh" \
    "$PDB_LIST" \
    "$MMCIF_DIR" \
    "$RESULT_DIR/structure_contacts.tsv" \
    "$RESULT_DIR/reordered.fasta" \
    "$RESULT_DIR" \
    "$RESULT_DIR/structure_sites.tsv" \
    "$COG_FASTA" \
    "$SITE_TITLE"
fi

DO_SUPERVISOR_SCAN="0"
if [[ -n "$PDB_LIST" && -f "$PDB_LIST" ]]; then
  if [[ "$RUN_SUPERVISOR_SCAN" == "1" ]]; then
    DO_SUPERVISOR_SCAN="1"
  elif [[ "$RUN_SUPERVISOR_SCAN" == "auto" && "$PDB_SOURCE" == "supervisor_mapping" ]]; then
    if command -v "$SUPERVISOR_PYTHON_BIN" >/dev/null 2>&1; then
      DO_SUPERVISOR_SCAN="1"
    fi
  fi
fi

if [[ "$DO_SUPERVISOR_SCAN" == "1" ]]; then
  SCAN_INPUT_DIR="$RESULT_DIR/pdb"
  if [[ ! -d "$SCAN_INPUT_DIR" ]]; then
    SCAN_INPUT_DIR="$RESULT_DIR/mmcif"
  fi
  bash "$PROJECT_ROOT/sh/06_supervisor_ligand_scan.sh" "$SCAN_INPUT_DIR" "$RESULT_DIR/supervisor_ligand_scan"
fi

echo "Done: $RESULT_DIR"
