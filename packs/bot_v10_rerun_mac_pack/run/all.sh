#!/usr/bin/env bash
# Full v10 rerun: the v8 270-cell scaling grid, then the v9 controls (A, B, C).
# From the pack folder:
#   bash run/all.sh
#   bash run/all.sh --smoke
# Re-run the same command to resume; finished cells are skipped.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PKG"

if [[ "${1:-}" != "" && "${1:-}" != "--smoke" ]]; then
  echo "usage: bash run/all.sh [--smoke]" >&2
  exit 2
fi

if [[ "$(uname)" == "Darwin" && -z "${LEX_AWAKE:-}" ]] && command -v caffeinate >/dev/null; then
  export LEX_AWAKE="caffeinate -dims (all.sh)"
  exec caffeinate -dims bash "$SCRIPT_DIR/all.sh" "$@"
fi
if [[ "$(uname)" == "Linux" && -z "${LEX_AWAKE:-}" ]]; then
  export LEX_AWAKE="none (Linux test; no caffeinate)"
fi
export LEX_AWAKE

if [[ "${1:-}" == "--smoke" ]]; then
  echo "SMOKE run: scaling (N = 2^16..2^18) then controls (N = 2^16 + small fake audio)."
else
  echo "v10 full rerun: scaling grid (~1 h on the M3 Max) then controls (~16 min)."
  echo "Expected runtime: about 1 hour 20 minutes."
fi
echo "Outputs: out/scaling/, out/controls/, out/env.txt, out/v10_mac_out.zip"
echo "Send back: out/v10_mac_out.zip"

# Each part finds Python and the compiler, builds if needed, and is resumable.
bash "$SCRIPT_DIR/scaling.sh" "$@"
bash "$SCRIPT_DIR/controls.sh" "$@"

# Combined zip of everything under out/, excluding the zip itself.
mkdir -p "$PKG/out"
rm -f "$PKG/out/v10_mac_out.zip"
# Also exclude the intermediate per-part zips so the combined archive stays small.
( cd "$PKG" && zip -qr out/v10_mac_out.zip out \
    -x out/v10_mac_out.zip out/scaling_mac_out.zip out/controls_mac_out.zip )
echo "Done. Send back: out/v10_mac_out.zip"
