#!/usr/bin/env bash
# Fill the paper's k grid: LEX (sig5) and wilber-interp at k=16,32,128,
# Ckmeans.1d.dp at k=32, plus a small-k check at k=4,8 on normal_n1048576.
# From the pack folder:
#   bash run/kgrid_fill.sh
#   bash run/kgrid_fill.sh --smoke
# Re-run the same command to resume; finished cells are skipped.
# Caps (optional, Ckmeans only): LEX_MEM_CAP_GIB=24 LEX_TIME_CAP_S=1800
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PKG"

SMOKE=0
for arg in "$@"; do
  case "$arg" in
    --smoke) SMOKE=1 ;;
    -h|--help)
      echo "usage: bash run/kgrid_fill.sh [--smoke]" >&2
      echo "  LEX_MEM_CAP_GIB   Ckmeans memory cap in GiB (default 24); cells over the cap are skipped" >&2
      echo "  LEX_TIME_CAP_S    Ckmeans wall-time cap in seconds (default 1800) per cell" >&2
      echo "  LEX_CKMEANS_METHOD  linear (default, fastest), loglinear, or quadratic" >&2
      exit 0
      ;;
    *)
      echo "usage: bash run/kgrid_fill.sh [--smoke]" >&2
      exit 2
      ;;
  esac
done

LEX_COMMIT=8a691ee005b657f985ffa299c57542b79570b08c
LEX_CPP_SHA=6c4a953c93914f73925ca501d2f29735c30d4c92aa14de82d98c8ced8db84a46

if [[ "$(uname)" == "Darwin" && -z "${LEX_AWAKE:-}" ]] && command -v caffeinate >/dev/null; then
  export LEX_AWAKE="caffeinate -dims (kgrid_fill.sh)"
  exec caffeinate -dims bash "$SCRIPT_DIR/kgrid_fill.sh" "$@"
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

mkdir -p "$PKG/data" "$PKG/work" "$PKG/out/lex" "$PKG/out/ckmeans" "$PKG/build" "$PKG/rlib"
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
    -I"$PKG/lex1d/cpp/include" -I"$PKG/lex1d/cpp/src" -I"$PKG/bench" -I"$PKG/bench/weighted_wilber" \
    "$PKG/bench/scale_cell.cpp" "$PKG/lex1d/cpp/src/"*.cpp \
    "$PKG/bench/weighted_wilber/interval_sum.cpp" \
    "$PKG/bench/weighted_wilber/kmeans_wilber.cpp" \
    -o "$BUILD/scale_cell"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -Wall -Wextra "$PKG/bench/f32_to_f64.cpp" -o "$BUILD/f32_to_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -Wall -Wextra "$PKG/bench/txt_to_f64.cpp" -o "$BUILD/txt_to_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra \
    "$PKG/bench/sort_f64.cpp" -o "$BUILD/sort_f64"
  $CXX_BIN -std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra \
    -I"$PKG/bench/weighted_wilber" \
    "$PKG/bench/score_sse.cpp" \
    "$PKG/bench/weighted_wilber/interval_sum.cpp" \
    "$PKG/bench/weighted_wilber/kmeans_wilber.cpp" \
    -o "$BUILD/score_sse"
fi

{
  echo "pack=bot_v1_kgrid_fill_pack"
  echo "smoke=$SMOKE"
  echo "lex1d_commit=$LEX_COMMIT"
  echo "lex_cpp_sha256=$LEX_CPP_SHA"
  echo "lex_method=sig5 (5-significant-digit grouping, k-preserving digit fallback, group-level Refine)"
  echo "exact_solver=wilber-interp (patched bench/weighted_wilber, interpolation search, unit weights)"
  echo "ckmeans_package=Ckmeans.1d.dp"
  echo "ckmeans_version=$CKMEANS_VER"
  echo "ckmeans_method=$CKMEANS_METHOD"
  echo "default_M=max(4096,32k)"
  echo "budget_note=M stays 4096 for k<=128 and grows with k above that; these cells evaluate the default budget policy"
  echo "k_lex_wilber=16,32,128 on uniform/normal/lognormal at N=2^20 and 2^23, and audio if present"
  echo "k_ckmeans=32 on the same datasets, sorted and raw"
  echo "k_small_check=4,8 on normal_n1048576 for LEX, wilber, and Ckmeans"
  echo "inputs=gen_synth.py SeedSequence(2048, spawn_key=(stream_id,N)); audio N=8168256 if present in data/"
  echo "protocol=1 warm-up + 3 timed reps, median (upper-middle); smoke uses 1 timed rep"
  echo "comparison=post-sort only; std::sort logged as sort_s and not added to post_s"
  echo "sse=two-pass long double raw_sse on the raw sorted values (scale_cell / score_sse)"
  echo "mem_cap_gib=$MEM_CAP"
  echo "time_cap_s=$TIME_CAP"
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
export LEX_MEM_CAP_GIB="${LEX_MEM_CAP_GIB:-24}"
export LEX_TIME_CAP_S="${LEX_TIME_CAP_S:-1800}"

if [[ "$SMOKE" == 1 ]]; then
  echo "SMOKE run: N=4096, three families; LEX and wilber at k in {16,32} (plus k=4,8 on normal); Ckmeans at k=32 (plus k=4,8 on normal); 1 warm-up + 1 timed rep."
else
  echo "k-grid fill: LEX (sig5) and wilber-interp at k=16,32,128; Ckmeans.1d.dp at k=32; small-k check at k=4,8 on normal_n1048576."
  echo "Expected runtime on an M3 Max: about 30 minutes (estimate), from earlier runs at nearby k and N."
  echo "Peak memory about 4.2 GiB (estimate) for Ckmeans k=32 at N≈8e6; LEX/wilber peak about 1.1 GiB (from earlier 8M wilber cells)."
  echo "The default 24 GiB Ckmeans cap is enough; you do not need to raise LEX_MEM_CAP_GIB."
fi
echo "Outputs: out/lex/, out/ckmeans/, out/env.txt, out/kgrid_fill_mac_out.zip"
echo "Send back: out/kgrid_fill_mac_out.zip"

py_args=()
if [[ "$SMOKE" == 1 ]]; then py_args+=(--smoke); fi

lex_rc=0
ck_rc=0
"$PYTHON_BIN" "$PKG/run/scaling_run.py" ${py_args[@]+"${py_args[@]}"} || lex_rc=$?
"$PYTHON_BIN" "$PKG/run/ckmeans_matched.py" ${py_args[@]+"${py_args[@]}"} || ck_rc=$?

zip_path="$PKG/out/kgrid_fill_mac_out.zip"
if [[ -f "$zip_path" ]]; then rm -f "$zip_path"; fi
(cd "$PKG" && zip -r out/kgrid_fill_mac_out.zip out/env.txt out/lex out/ckmeans \
  -x "out/kgrid_fill_mac_out.zip" -x "out/*/*.zip")

if [[ "$lex_rc" -ne 0 || "$ck_rc" -ne 0 ]]; then
  echo "error: fill finished with lex_rc=$lex_rc ck_rc=$ck_rc. Partial results are in out/ and out/kgrid_fill_mac_out.zip." >&2
  exit 1
fi
echo "Done. Results: out/lex/summary.md, out/ckmeans/summary.md, and out/kgrid_fill_mac_out.zip"
