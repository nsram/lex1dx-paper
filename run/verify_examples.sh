#!/usr/bin/env bash
# Recheck the exact-fraction worked examples in the paper.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/verify_examples.py
