#!/usr/bin/env bash
# Redraw paper/figs/v5_seven.pdf (r6 Figure 2) from the exact-fraction example.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/v5_seven.py
