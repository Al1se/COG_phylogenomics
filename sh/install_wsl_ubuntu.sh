#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v apt >/dev/null 2>&1; then
  echo "This installer is intended for Ubuntu/Debian under WSL."
  echo "Install dependencies manually on your system."
  exit 1
fi

sudo apt update
sudo apt install -y \
  python3 \
  python3-venv \
  python3-pip \
  muscle \
  iqtree \
  git \
  curl

bash "$PROJECT_ROOT/sh/install_python_env.sh"

cat <<EOF
WSL/Ubuntu setup completed.

Next steps:
  source "$PROJECT_ROOT/.venv/bin/activate"
  ls "$PROJECT_ROOT/input"

Expected input files:
  cog-24.cog.csv
  cog-24.org.csv
  cog-24.tax.csv
  2296Genomes.prot.fasta

First test run:
  bash "$PROJECT_ROOT/sh/run_full_analysis.sh" COG0086 "$PROJECT_ROOT/results/COG0086" short
EOF
