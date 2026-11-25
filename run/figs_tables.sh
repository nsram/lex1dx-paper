#!/usr/bin/env bash
# Regenerate paper/figs/v4_*.pdf and paper/tables/v4_*.tex from the packaged CSVs.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/v4_figs_tables.py
