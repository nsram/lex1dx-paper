#!/usr/bin/env bash
# Rerun LEX sig5 k-grid cells on the reported path (Refine on precision-group
# atoms; map to raw indices only for SSE). Does not rerun wilber or Ckmeans.
# From the pack folder:
#   bash run/kgrid_lex_fix.sh
#   bash run/kgrid_lex_fix.sh --smoke
# Re-run the same command to resume; finished cells are skipped.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PKG"

SMOKE=0
for arg in "$@"; do
  case "$arg" in
    --smoke) SMOKE=1 ;;
    -h|--help)
      echo "usage: bash run/kgrid_lex_fix.sh [--smoke]" >&2
      exit 0
      ;;
    *)
      echo "usage: bash run/kgrid_lex_fix.sh [--smoke]" >&2
      exit 2
      ;;
  esac
done

LEX_COMMIT=8a691ee005b657f985ffa299c57542b79570b08c
LEX_CPP_SHA=6c4a953c93914f73925ca501d2f29735c30d4c92aa14de82d98c8ced8db84a46
# Reported LEX harness: bot_v8_scaling_mac_pack/bench/scale_cell.cpp, unmodified.
SCALE_CELL_SHA=3b5f0e910f486d9968283483b0ba7b61b4e48b002e31651abf46ae10010a80e6
TUNED_COLLAPSE_SHA=98839e35ef82c40672a0e70cb15f910dcddca024dd4b4575aa36f7bc130d3810

if [[ "$(uname)" == "Darwin" && -z "${LEX_AWAKE:-}" ]] && command -v caffeinate >/dev/null; then
  export LEX_AWAKE="caffeinate -dims (kgrid_lex_fix.sh)"
  exec caffeinate -dims bash "$SCRIPT_DIR/kgrid_lex_fix.sh" "$@"
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

LEX_CPP="$PKG/lex1d/cpp/src/lex.cpp"
if [[ ! -f "$LEX_CPP" ]]; then
  echo "error: missing lex1d sources at lex1d/cpp/src/lex.cpp" >&2
  exit 2
fi
got=$(sha256_of "$LEX_CPP")
if [[ "$got" != "$LEX_CPP_SHA" ]]; then
  echo "error: lex1d/cpp/src/lex.cpp SHA-256 is $got, expected $LEX_CPP_SHA (commit $LEX_COMMIT)" >&2
  exit 2
fi

SCALE_CPP="$PKG/bench/scale_cell.cpp"
got=$(sha256_of "$SCALE_CPP")
if [[ "$got" != "$SCALE_CELL_SHA" ]]; then
  echo "error: bench/scale_cell.cpp SHA-256 is $got, expected $SCALE_CELL_SHA (reported LEX harness, unmodified)" >&2
  exit 2
fi
COLLAPSE="$PKG/bench/tuned_collapse.hpp"
got=$(sha256_of "$COLLAPSE")
if [[ "$got" != "$TUNED_COLLAPSE_SHA" ]]; then
  echo "error: bench/tuned_collapse.hpp SHA-256 is $got, expected $TUNED_COLLAPSE_SHA" >&2
  exit 2
fi

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

mkdir -p "$PKG/data" "$PKG/work" "$PKG/out/lex" "$PKG/build"

BUILD="$PKG/build"
if command -v cmake >/dev/null; then
  cmake -S "$PKG" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$CXX_BIN"
  cmake --build "$BUILD" --config Release
else
  echo "cmake not found; compiling directly with $CXX_BIN" >&2
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra \
    -I"$PKG/lex1d/cpp/include" -I"$PKG/lex1d/cpp/src" -I"$PKG/bench" -I"$PKG/bench/weighted_wilber" \
    "$PKG/bench/scale_cell.cpp" "$PKG/lex1d/cpp/src/"*.cpp \
    "$PKG/bench/weighted_wilber/interval_sum.cpp" \
    "$PKG/bench/weighted_wilber/kmeans_wilber.cpp" \
    -o "$BUILD/scale_cell"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -Wall -Wextra "$PKG/bench/f32_to_f64.cpp" -o "$BUILD/f32_to_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -Wall -Wextra "$PKG/bench/txt_to_f64.cpp" -o "$BUILD/txt_to_f64"
fi

{
  echo "pack=bot_v2_kgrid_fill_pack"
  echo "smoke=$SMOKE"
  echo "lex1d_commit=$LEX_COMMIT"
  echo "lex_cpp_sha256=$LEX_CPP_SHA"
  echo "scale_cell_sha256=$SCALE_CELL_SHA"
  echo "tuned_collapse_sha256=$TUNED_COLLAPSE_SHA"
  echo "lex_method=sig5 (always 5-significant-digit grouping; Refine on precision-group atoms; map to raw indices only for SSE)"
  echo "default_M=max(4096,32k)"
  echo "budget_note=M stays 4096 for k<=128 and grows with k above that; these cells evaluate the default budget policy"
  echo "k_lex=16,32,128 on uniform/normal/lognormal at N=2^20 and 2^23, and audio if present"
  echo "k_small_check=4,8 on normal_n1048576"
  echo "not_rerun=wilber-interp, Ckmeans.1d.dp"
  echo "inputs=gen_synth.py SeedSequence(2048, spawn_key=(stream_id,N)); audio N=8168256 if present in data/"
  echo "protocol=1 warm-up + 3 timed reps, median (upper-middle); smoke uses 1 timed rep"
  echo "comparison=post-sort only; std::sort logged as sort_s and not added to post_s"
  echo "sse=two-pass long double raw_sse on the raw sorted values (scale_cell)"
  echo "date=$(date)"
  echo "uname=$(uname -a)"
  if [[ "$(uname)" == "Darwin" ]] && command -v sw_vers >/dev/null; then sw_vers; fi
  echo "python=$PYTHON_BIN"
  "$PYTHON_BIN" --version 2>&1 | sed 's/^/python_version=/'
  "$PYTHON_BIN" -c 'import numpy; print("numpy="+numpy.__version__)' 2>/dev/null || true
  echo "cxx=$CXX_BIN"
  "$CXX_BIN" --version 2>&1 | sed 's/^/cxx_version=/'
  if command -v cmake >/dev/null; then cmake --version 2>&1 | head -1 | sed 's/^/cmake=/'; else echo "cmake=(none, direct compile)"; fi
  echo "awake=${LEX_AWAKE:-}"
  echo "omp_num_threads=${OMP_NUM_THREADS}"
  echo "cwd=$PKG"
} > "$PKG/out/env.txt"

export PYTHON_BIN

if [[ "$SMOKE" == 1 ]]; then
  echo "SMOKE run: N=4096, three families; LEX sig5 at k in {16,32} (plus k=4,8 on normal); 1 warm-up + 1 timed rep."
else
  echo "LEX sig5 k-grid fix: k=16,32,128 on the 1M and 8M synthetics and audio, plus k=4,8 on normal_n1048576."
  echo "Expected runtime on an M3 Max: about 5 minutes (estimate), including first-time data generation."
  echo "Peak memory about 130 MiB (estimate), from earlier 8M sig5 cells on this code path."
fi
echo "Outputs: out/lex/, out/env.txt, out/kgrid_lex_fix_mac_out.zip"
echo "Send back: out/kgrid_lex_fix_mac_out.zip"

py_args=()
if [[ "$SMOKE" == 1 ]]; then py_args+=(--smoke); fi
"$PYTHON_BIN" "$PKG/run/scaling_run.py" ${py_args[@]+"${py_args[@]}"}

zip_path="$PKG/out/kgrid_lex_fix_mac_out.zip"
if [[ -f "$zip_path" ]]; then rm -f "$zip_path"; fi
(cd "$PKG" && zip -r out/kgrid_lex_fix_mac_out.zip out/env.txt out/lex \
  -x "out/kgrid_lex_fix_mac_out.zip" -x "out/*/*.zip")
echo "Done. Results: out/lex/summary.md and out/kgrid_lex_fix_mac_out.zip"
