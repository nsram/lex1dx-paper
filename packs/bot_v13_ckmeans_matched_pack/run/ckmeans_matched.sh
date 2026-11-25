#!/usr/bin/env bash
# Matched Ckmeans.1d.dp rerun on the lex1dx v10 inputs and protocol.
# From the pack folder:
#   bash run/ckmeans_matched.sh
#   bash run/ckmeans_matched.sh --smoke
#   bash run/ckmeans_matched.sh --sorted-only
# Re-run the same command to resume; finished cells are skipped.
# Caps (optional): LEX_MEM_CAP_GIB=24 LEX_TIME_CAP_S=1800
# Method (optional): LEX_CKMEANS_METHOD=linear   # default; fastest exact method
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PKG"

SMOKE=0
SORTED_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --smoke) SMOKE=1 ;;
    --sorted-only) SORTED_ONLY=1 ;;
    -h|--help)
      echo "usage: bash run/ckmeans_matched.sh [--smoke] [--sorted-only]" >&2
      echo "  LEX_MEM_CAP_GIB   memory cap in GiB (default 24); cells over the cap are skipped" >&2
      echo "  LEX_TIME_CAP_S    wall-time cap in seconds (default 1800) per cell" >&2
      echo "  LEX_CKMEANS_METHOD  linear (default, fastest), loglinear, or quadratic" >&2
      exit 0
      ;;
    *)
      echo "usage: bash run/ckmeans_matched.sh [--smoke] [--sorted-only]" >&2
      exit 2
      ;;
  esac
done

if [[ "$(uname)" == "Darwin" && -z "${LEX_AWAKE:-}" ]] && command -v caffeinate >/dev/null; then
  export LEX_AWAKE="caffeinate -dims (ckmeans_matched.sh)"
  exec caffeinate -dims bash "$SCRIPT_DIR/ckmeans_matched.sh" "$@"
fi
if [[ "$(uname)" == "Linux" && -z "${LEX_AWAKE:-}" ]]; then
  export LEX_AWAKE="none (Linux test; no caffeinate)"
fi

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

sha256_of() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

find_python() {
  if [[ -n "${PYTHON:-}" ]] && command -v "$PYTHON" >/dev/null; then
    command -v "$PYTHON"
    return 0
  fi
  local p
  for p in python3 python3.13 python3.12 python3.11 python3.10 python; do
    if command -v "$p" >/dev/null; then
      command -v "$p"
      return 0
    fi
  done
  return 1
}

find_cxx() {
  if [[ -n "${CXX:-}" ]] && command -v "$CXX" >/dev/null; then
    command -v "$CXX"
    return 0
  fi
  local c
  for c in c++ clang++ g++ /usr/bin/c++ /usr/bin/clang++; do
    if command -v "$c" >/dev/null; then
      command -v "$c"
      return 0
    fi
  done
  return 1
}

find_rscript() {
  if [[ -n "${RSCRIPT:-}" ]] && command -v "$RSCRIPT" >/dev/null; then
    command -v "$RSCRIPT"
    return 0
  fi
  local p
  for p in Rscript \
           /usr/local/bin/Rscript \
           /opt/homebrew/bin/Rscript \
           /Library/Frameworks/R.framework/Resources/bin/Rscript \
           /usr/bin/Rscript; do
    if command -v "$p" >/dev/null; then
      command -v "$p"
      return 0
    fi
  done
  return 1
}

if ! PYTHON_BIN=$(find_python); then
  echo "error: no Python found (tried python3, python). Install Python 3.8+ and NumPy 1.17+." >&2
  exit 2
fi
if ! "$PYTHON_BIN" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 8) else 2)'; then
  echo "error: $($PYTHON_BIN --version 2>&1) is too old; need Python 3.8+." >&2
  exit 2
fi
if ! "$PYTHON_BIN" -c 'import numpy' 2>/dev/null; then
  echo "error: NumPy 1.17+ (PCG64 / default_rng) is required to generate the synthetic data." >&2
  echo "  $PYTHON_BIN -m pip install 'numpy>=1.17'" >&2
  exit 2
fi
if ! CXX_BIN=$(find_cxx); then
  echo "error: no C++ compiler found (tried c++, clang++, g++). On a Mac, install the Xcode command line tools." >&2
  exit 2
fi
if ! RSCRIPT_BIN=$(find_rscript); then
  echo "error: no Rscript found. Install R (https://cran.r-project.org) and make sure Rscript is on PATH." >&2
  exit 2
fi

mkdir -p "$PKG/data" "$PKG/work" "$PKG/out/ckmeans_matched" "$PKG/build" "$PKG/rlib"
export R_LIBS="$PKG/rlib${R_LIBS:+:$R_LIBS}"
export R_LIBS_USER="$PKG/rlib"

# Install Ckmeans.1d.dp into the pack-local library if it is not already there.
if ! "$RSCRIPT_BIN" -e 'if (!requireNamespace("Ckmeans.1d.dp", quietly=TRUE)) quit(status=2)' >/dev/null 2>&1; then
  echo "Installing Ckmeans.1d.dp from CRAN into rlib/ ..."
  "$RSCRIPT_BIN" -e 'install.packages("Ckmeans.1d.dp", lib="'"$PKG/rlib"'", repos="https://cloud.r-project.org", dependencies=TRUE)'
  if ! "$RSCRIPT_BIN" -e 'if (!requireNamespace("Ckmeans.1d.dp", quietly=TRUE)) quit(status=2)' >/dev/null 2>&1; then
    echo "error: failed to install Ckmeans.1d.dp. On a Mac, install the Xcode command line tools and try again." >&2
    exit 2
  fi
fi

CKMEANS_VER=$("$RSCRIPT_BIN" -e 'cat(as.character(packageVersion("Ckmeans.1d.dp")))')
CKMEANS_METHOD="${LEX_CKMEANS_METHOD:-linear}"
MEM_CAP="${LEX_MEM_CAP_GIB:-24}"
TIME_CAP="${LEX_TIME_CAP_S:-1800}"

BUILD="$PKG/build"
if command -v cmake >/dev/null; then
  cmake -S "$PKG" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$CXX_BIN"
  cmake --build "$BUILD" --config Release
else
  echo "cmake not found; compiling directly with $CXX_BIN" >&2
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra \
    "$PKG/bench/sort_f64.cpp" -o "$BUILD/sort_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -Wall -Wextra \
    "$PKG/bench/txt_to_f64.cpp" -o "$BUILD/txt_to_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra \
    -I"$PKG/bench/weighted_wilber" \
    "$PKG/bench/score_sse.cpp" \
    "$PKG/bench/weighted_wilber/interval_sum.cpp" \
    "$PKG/bench/weighted_wilber/kmeans_wilber.cpp" \
    -o "$BUILD/score_sse"
fi

{
  echo "pack=bot_v11_ckmeans_matched_pack"
  echo "smoke=$SMOKE"
  echo "sorted_only=$SORTED_ONLY"
  echo "ckmeans_package=Ckmeans.1d.dp"
  echo "ckmeans_version=$CKMEANS_VER"
  echo "ckmeans_method=$CKMEANS_METHOD"
  echo "mem_cap_gib=$MEM_CAP"
  echo "time_cap_s=$TIME_CAP"
  echo "protocol=1 warm-up + 3 timed reps, median (upper-middle); smoke uses 1 timed rep"
  echo "sse=v10 scale_cell raw_sse (two-pass long double on raw sorted values)"
  echo "exact_reference=wilber-interp (patched bench/weighted_wilber, interpolation search)"
  echo "sort=std::sort with x=x+0.0, timed separately as sort_s; Ckmeans internal sort is not excluded"
  echo "inputs=v10 gen_synth.py SeedSequence(2048, spawn_key=(stream_id,N)); audio N=8168256"
  echo "cells_full=uniform,normal,lognormal at N=1048576 and N=8388608, k in {16,128,512}; audio if present"
  echo "date=$(date)"
  echo "uname=$(uname -a)"
  if [[ "$(uname)" == "Darwin" ]] && command -v sw_vers >/dev/null; then sw_vers; fi
  echo "python=$PYTHON_BIN"
  "$PYTHON_BIN" --version 2>&1 | sed 's/^/python_version=/'
  "$PYTHON_BIN" -c 'import numpy; print("numpy="+numpy.__version__)' 2>/dev/null || true
  echo "cxx=$CXX_BIN"
  "$CXX_BIN" --version 2>&1 | sed 's/^/cxx_version=/'
  echo "rscript=$RSCRIPT_BIN"
  "$RSCRIPT_BIN" --version 2>&1 | sed 's/^/r_version=/'
  if command -v cmake >/dev/null; then cmake --version 2>&1 | head -1 | sed 's/^/cmake=/'; else echo "cmake=(none, direct compile)"; fi
  echo "awake=${LEX_AWAKE:-}"
  echo "omp_num_threads=${OMP_NUM_THREADS}"
  echo "cwd=$PKG"
} > "$PKG/out/env.txt"

export PYTHON_BIN
export RSCRIPT_BIN
export LEX_CKMEANS_METHOD="$CKMEANS_METHOD"
# Caps already in the environment if the user set them; export defaults so Python sees them.
export LEX_MEM_CAP_GIB="${LEX_MEM_CAP_GIB:-24}"
export LEX_TIME_CAP_S="${LEX_TIME_CAP_S:-1800}"

if [[ "$SMOKE" == 1 ]]; then
  echo "SMOKE run: N=4096, three families, k in {4,16}, 1 warm-up + 1 timed rep."
else
  echo "Ckmeans.1d.dp matched rerun (method=$CKMEANS_METHOD, version=$CKMEANS_VER)."
  echo "Expected runtime on the M3 Max: about 45 minutes (both sorted and unsorted fits;"
  echo "k=512 at N≈8e6 is skipped under the default 24 GiB memory cap)."
  echo "About 25 minutes with --sorted-only."
fi
echo "Outputs: out/ckmeans_matched/, out/env.txt, out/ckmeans_matched_mac_out.zip"
echo "Send back: out/ckmeans_matched_mac_out.zip"

py_args=()
if [[ "$SMOKE" == 1 ]]; then py_args+=(--smoke); fi
if [[ "$SORTED_ONLY" == 1 ]]; then py_args+=(--sorted-only); fi
"$PYTHON_BIN" "$PKG/run/ckmeans_matched.py" ${py_args[@]+"${py_args[@]}"}
