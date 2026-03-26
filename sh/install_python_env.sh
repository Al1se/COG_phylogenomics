#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${1:-.venv}"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$PYTHON_BIN" -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel
"$VENV_DIR/bin/pip" install -r "$PROJECT_ROOT/requirements.txt"

cat <<EOF
Python environment created at: $VENV_DIR

Activate it with:
  source "$VENV_DIR/bin/activate"

Installed Python dependencies from:
  $PROJECT_ROOT/requirements.txt

External tools are still required separately:
  - MUSCLE: required for MSA steps
  - IQ-TREE: required for tree building in the main pipeline
  - PyMOL: optional, only for manual structure inspection
  - Python 2: optional, only for supervisor legacy tools
EOF
