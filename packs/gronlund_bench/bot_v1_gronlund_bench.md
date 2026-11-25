# Grønlund et al. authors' code vs LEX — small bench (bot_v1)

Work directory: `/workspace/lexicon/work/gronlund_bench`. Nothing under `/workspace/lexicon/originals` was edited. The cloned `kmeans1d` sources were not patched.

Run window: 2026-11-24 07:57:50–08:17:12 America/New_York. Total wall time of `run_bench.py` was 1142.2 s (correctness plus scale). Scale phase 915.6 s; no scale cell was skipped.

## Environment

```
OS: Debian GNU/Linux 13 (trixie), Linux 6.12.113+ x86_64
CPU: Intel(R) Xeon(R) Processor, 8 cores, 2400.000 MHz (from /proc/cpuinfo)
      hypervisor KVM; L1d 384 KiB, L1i 256 KiB, L2 16 MiB, L3 320 MiB
RAM: 15 GiB (no swap)
compiler: g++ (Debian 14.2.0-19) 14.2.0
cmake: 3.31.6
make: GNU Make 4.4.1
GNU time: /usr/bin/time (package time 1.9-0.2)
python: 3.13 / numpy 2.5.5
```

This box is not the M3 Max used for the LEX paper numbers. Only comparisons within this box are on the same machine. `Ckmeans.1d.dp` (R) was not installed. The exact reference is therefore LEX `lex_cli --preset exact` (`lex::ckmeans_exact`).

## Repo identity

```
git clone https://github.com/gronlund/kmeans1d.git
commit: 585386d4678000691706acdb954911be8fea1c5c
date:   2021-04-09T12:39:29+02:00
        = 2021-04-09 06:39:29 America/New_York
author: Allan Grønlund <allan.g.joergensen@gmail.com>
subject: small bug fixes
LICENSE: MIT (LICENSE file)
README: "The code is used for the experiments in the paper
        Fast Exact k-Means, k-Medians and Bregman Divergence Clustering in 1D
        (https://arxiv.org/abs/1701.07204)"
```

## Algorithms in this repository

From `kmeans.h` and the `.cpp` files (class `name()` in parentheses):

| class | `name()` | What the source implements | Space in `compute()` | Path / centers |
|---|---|---|---|---|
| `kmeans_linear` | `linear` | SMAWK / "matrix searching" 1-D k-means; comment in `kmeans_linear.cpp`: "This code implements the matrix searching algorithm for computing k-means." O(kn) time fill. | two rows (`row`, `row_prev`), O(n) | `compute()` returns cost only (`row[n-1]`). `compute_and_report()` calls `kmeans_dp::report`, a Hirschberg-style split (solve left `k/2` and reversed right `k-k/2`, recurse) and then `compute()` again for the cost. |
| `kmeans_monotone` | `monotone` | Monotone / quadrangle-inequality divide-and-conquer row fill (`fill_row_rec`). | two rows, O(n) | same `report` as above |
| `kmeans_slow` | `slow` | Naive DP fill, O(kn²). | two rows, O(n) | same `report` |
| `kmeans_hirschberg_larmore` | `hirc` | Hirschberg–Larmore least-weight subsequence (`basic`) plus a search on a penalty `lambda` until the path has k segments. O(n) arrays `f`, `bestleft`. This is a penalty / Lagrangian method, not the 1975 Hirschberg DP reconstruction. | O(n) | `compute()` already reconstructs centers via `get_actual_cost`. `compute_and_report()` is the same call. |
| `kmeans_wilber` | `wilber` | Wilber / SMAWK unconstrained (penalized) problem, then search on `lambda`. Two strategies: `BINARY` and `INTERPOLATION` (the default in `timing.cpp`). | O(n) | `compute()` reconstructs centers and `path` in `get_actual_cost`. |
| `kmeans_lloyd_slow`, `kmeans_lloyd_fast` | `lloyd` | Lloyd heuristic. Not run (not exact). | — | — |

There is no implementation in this repo that stores a full k×n back-pointer table. The O(kn)-time DPs are `linear` (SMAWK) and `monotone` (divide-and-conquer). The O(n)-space / table-free variants that return a clustering are: `linear`/`monotone`/`slow` with `compute_and_report` (Hirschberg-style recomputation), `hirc`, `wilber-binary`, and `wilber-interp`.

The library assumes the input is already sorted. The driver sorts a copy and records whether the file was sorted.

`cost` is the L2 interval SSE from `interval_sum::cost_interval_l2` (sum of squared deviations from the segment mean).

## Build

`make` and `cmake` were not initially installed; they were installed with `sudo apt-get install make cmake time` (GNU `time` was also missing). No source changes were required after that.

### Grønlund (`make timing`), no local patch

```
cd /workspace/lexicon/work/gronlund_bench/kmeans1d
make timing
```

Compiler and flags (from the Makefile, `T ?= RELEASE`):

```
g++ -Wall -Wextra -fPIE -O2 -std=c++11 -c -o <file>.o <file>.cpp
g++ -fopenmp -Wall -Wextra -fPIE -O2 -std=c++11 -o timing <objs> timing.cpp
```

Exit status 0. The binary `timing` was produced. There were no errors. There were 17 warnings (unused variables/parameters, a signed/unsigned comparison in `kmeans_monotone.cpp`, and in `timing.cpp` the string-literal comparisons `argv[1] == "uniform"`). Full log: `build_timing.log`. Asserts are left enabled (`-O2`, no `-DNDEBUG`).

The stock `timing` binary generates its own inputs and times several algorithms; it was compiled to confirm the stock build, but the measurements below use a separate driver linked against the same object files.

### Driver (written for this bench, not part of the upstream repo)

```
g++ -Wall -Wextra -O2 -std=c++11 -I kmeans1d -o gronlund_driver driver.cpp \
  kmeans1d/kmeans_dp.o kmeans1d/kmeans_slow.o kmeans1d/kmeans_linear.o \
  kmeans1d/kmeans_monotone.o kmeans1d/kmeans_lloyd.o \
  kmeans1d/kmeans_hirschberg_larmore.o kmeans1d/kmeans_wilber.o \
  kmeans1d/interval_sum.o
```

Usage: `gronlund_driver --alg ALG --k K [--mode compute|report] [--out FILE] DATAFILE`.

`ALG` is one of `linear`, `monotone`, `slow`, `hirc`, `wilber-binary`, `wilber-interp`. `mode=compute` calls `compute()`; `mode=report` calls `compute_and_report()`. Wall time inside the driver (`driver_s`) is `std::chrono::steady_clock` around that call only, after I/O and after sorting.

### LEX (paper build)

LEX was built out of tree from the paper freeze, not from a later tree:

- source: `/workspace/lexicon/work/package_v25/bot_v25_lex_pack/lex1d` (the pack states this is commit `8a691ee005b657f985ffa299c57542b79570b08c`, verbatim)
- `lex.cpp` SHA-256: `6c4a953c93914f73925ca501d2f29735c30d4c92aa14de82d98c8ced8db84a46` (matches the pack README)
- build dir: `/workspace/lexicon/work/gronlund_bench/lex_build`

```
cmake -S /workspace/lexicon/work/package_v25/bot_v25_lex_pack/lex1d/cpp \
      -B /workspace/lexicon/work/gronlund_bench/lex_build \
      -DCMAKE_BUILD_TYPE=Release -DLEX_BUILD_TESTS=OFF -DLEX_BUILD_CLI=ON
cmake --build /workspace/lexicon/work/gronlund_bench/lex_build -j8
```

Flags recorded in `lex_build/CMakeFiles/lex.dir/flags.make`: `-O3 -DNDEBUG -std=c++17 -fPIC -ffp-contract=off -Wall -Wextra -Wpedantic`.

Scale runs used the paper default `M = 4096`:

```
lex_cli FILE K --preset hybrid --coarse 4096
```

Correctness used the exact path:

```
lex_cli FILE K --preset exact
```

## Data

Generated by `gen_data.py` into `data/`. NumPy PCG64, one double per line, `%.17g`.

| file | N | unique values | generator |
|---|---:|---:|---|
| `uniform_n{1000,10000}.txt` | 1e3 / 1e4 | = N | `Generator(PCG64(1)).uniform(0, 1)` |
| `normal_n{1000,10000,1e5,1e6}.txt` | those N | = N | `Generator(PCG64(2)).normal(0, 1)` |
| `lognormal_n{1000,10000}.txt` | 1e3 / 1e4 | = N | `Generator(PCG64(3)).lognormal(0, 1)` |
| `ties_n{1000,10000}.txt` | 1e3 / 1e4 | 20 | `Generator(PCG64(4)).integers(0, 20)` as float (unsorted) |
| `audio_prefix_n{1e5,1e6}.txt` | 1e5 / 1e6 | 99,771 / 989,460 | first N lines of the paper audio file `dharini-abhogi.txt` (8,168,256 rows). There was no smaller pre-made audio copy. This is a prefix, not a paper-defined slice. Source read only: `/workspace/lexicon/originals/lex_zip/extracted/lex/empirical/input/dharini-abhogi.txt`. |

For N = 1e3, D ≤ 1000 ≤ M = 4096, so `lex_hybrid` would also have been exact; the reference used is still `--preset exact`. For N = 1e4, D = 10,000 > 4096 except on the ties files (D = 20), so the exact preset is required. On the scale sets, D > 4096, so LEX hybrid is not exact.

## How each cell was run

```
/usr/bin/time -v -o LOG.time timeout --signal=KILL T \
  ./gronlund_driver --alg ALG --k K --mode MODE --out LOG.result DATA
# or
/usr/bin/time -v -o LOG.time timeout --signal=KILL T \
  ./lex_build/lex_cli DATA K --preset exact|hybrid [--coarse 4096]
```

Correctness timeout T = 60 s; scale T = 180 s. `elapsed s` and peak RSS are from GNU `/usr/bin/time -v` (elapsed includes process start, file parse, and, for the Grønlund driver, the sort). Peak RSS is `Maximum resident set size (kbytes) / 1024` in MiB. `driver s` is the in-process timer around `compute` / `compute_and_report` only. LEX has no `driver s`. Orchestration: `python3 -u run_bench.py`. Logs: `logs/`. CSVs: `results/correctness.csv`, `results/scale.csv`.

A smoke test on `{1,2,3,4,5,6,7,8}`, k = 2, all six exact Grønlund algorithms and LEX exact, all returned SSE = 10 and centers 2.5, 6.5.

## Correctness (N ∈ {1e3, 1e4}, k ∈ {4, 16, 64})

Relative difference is `(SSE_method − SSE_exact) / SSE_exact`, where `SSE_exact` is `lex_cli --preset exact`. All methods were run in `report` mode so centers were requested.

Table: relative SSE vs LEX exact. ABORT = SIGABRT (assert). TIMEOUT = killed at 60 s (SIGKILL, `timeout` exit 137).

| dataset | N | k | linear | monotone | slow | hirc | wilber-binary | wilber-interp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| uniform_n1000 | 1000 | 4 | 4.11e-14 | 4.11e-14 | 4.11e-14 | 4.11e-14 | 4.11e-14 | 4.11e-14 |
| uniform_n1000 | 1000 | 16 | 3.43e-13 | 3.43e-13 | 3.43e-13 | 3.43e-13 | 3.43e-13 | 3.43e-13 |
| uniform_n1000 | 1000 | 64 | 6.97e-12 | 6.97e-12 | 6.97e-12 | 6.97e-12 | 6.97e-12 | 6.97e-12 |
| normal_n1000 | 1000 | 4 | 3.27e-15 | 3.27e-15 | 3.27e-15 | 3.27e-15 | 3.27e-15 | 3.27e-15 |
| normal_n1000 | 1000 | 16 | 7.18e-14 | 7.18e-14 | 7.18e-14 | 7.18e-14 | 7.18e-14 | 7.18e-14 |
| normal_n1000 | 1000 | 64 | 1.40e-12 | 1.40e-12 | 1.40e-12 | 1.40e-12 | 1.40e-12 | 1.40e-12 |
| lognormal_n1000 | 1000 | 4 | 1.88e-14 | 1.88e-14 | 1.88e-14 | 1.88e-14 | 1.88e-14 | 1.88e-14 |
| lognormal_n1000 | 1000 | 16 | 3.76e-13 | 3.76e-13 | 3.76e-13 | 3.76e-13 | 3.76e-13 | 3.76e-13 |
| lognormal_n1000 | 1000 | 64 | 1.12e-11 | 1.12e-11 | 1.12e-11 | 1.12e-11 | 1.12e-11 | 1.12e-11 |
| ties_n1000 | 1000 | 4 | 6.34e-15 | 6.34e-15 | 6.34e-15 | 6.34e-15 | 6.34e-15 | 6.34e-15 |
| ties_n1000 | 1000 | 16 | -1.12e-15 | -1.12e-15 | -1.12e-15 | -1.12e-15 | -1.12e-15 | -1.12e-15 |
| ties_n1000 | 1000 | 64 | -1.0† | -1.0† | -1.0† | ABORT | TIMEOUT | -1.0† |
| uniform_n10000 | 10000 | 4 | -1.48e-13 | -1.48e-13 | -1.48e-13 | -1.48e-13 | -1.48e-13 | -1.48e-13 |
| uniform_n10000 | 10000 | 16 | -2.72e-12 | -2.72e-12 | -2.72e-12 | -2.72e-12 | -2.72e-12 | -2.72e-12 |
| uniform_n10000 | 10000 | 64 | -4.58e-11 | -4.58e-11 | -4.58e-11 | -4.58e-11 | -4.58e-11 | -4.58e-11 |
| normal_n10000 | 10000 | 4 | -5.22e-14 | -5.22e-14 | -5.22e-14 | -5.22e-14 | -5.22e-14 | -5.22e-14 |
| normal_n10000 | 10000 | 16 | -7.21e-13 | -7.21e-13 | -7.21e-13 | -7.21e-13 | -7.21e-13 | -7.21e-13 |
| normal_n10000 | 10000 | 64 | -1.09e-11 | -1.09e-11 | -1.09e-11 | -1.09e-11 | -1.09e-11 | -1.09e-11 |
| lognormal_n10000 | 10000 | 4 | -1.44e-14 | -1.44e-14 | -1.44e-14 | -1.44e-14 | -1.44e-14 | -1.44e-14 |
| lognormal_n10000 | 10000 | 16 | -2.36e-13 | -2.36e-13 | -2.36e-13 | -2.36e-13 | -2.36e-13 | -2.36e-13 |
| lognormal_n10000 | 10000 | 64 | -7.06e-12 | -7.06e-12 | -7.06e-12 | -7.06e-12 | -7.06e-12 | -7.06e-12 |
| ties_n10000 | 10000 | 4 | 4.85e-14 | 4.85e-14 | 4.85e-14 | 4.85e-14 | 4.85e-14 | 4.85e-14 |
| ties_n10000 | 10000 | 16 | 3.24e-14 | 3.24e-14 | 3.24e-14 | 3.24e-14 | 3.24e-14 | 3.24e-14 |
| ties_n10000 | 10000 | 64 | -1.0† | -1.0† | -1.0† | ABORT | TIMEOUT | ABORT |

† On `ties_*` k = 64 there are 20 distinct values and k > D. LEX exact reports 20 clusters (`atoms=20`) and a residual SSE of 1.6362552896094227e-26 (N = 1000) or 1.1696982983584281e-26 (N = 10,000). The successful Grønlund methods report SSE = 0 exactly, so `(0 − ε)/ε = −1`. That is a comparison of a numerical zero to an exact zero, not a large clustering error.

On those same two cells, `linear` / `monotone` / `slow` `report` returned `n_centers = 23` (not 64 and not 20). `wilber-interp` on `ties_n1000` k = 64 returned `n_centers = 64` and SSE = 0. On every other successful correctness cell (k ≤ D), `n_centers = k`.

Among the 132 successful comparisons with LEX exact SSE > 1e-20, the maximum absolute relative SSE difference was **4.578706987184621 × 10⁻¹¹**, on `uniform_n10000` k = 64 (all six methods; SSE 0.194028678129202 vs exact 0.194028678138086). The largest positive relative difference in that set was 1.1215143007117075 × 10⁻¹¹ (`lognormal_n1000` k = 64). On each successful cell the six Grønlund methods reported the same SSE as each other to the printed precision.

## Crashes and timeouts (correctness only)

All five failures are on the tie-heavy integer sets with k = 64 > D = 20. No scale cell failed.

1. **`hirc` on `ties_n1000` k = 64 and `ties_n10000` k = 64** — SIGABRT, rc = 134, GNU time "Command terminated by signal 6", elapsed 0.01 s.

   ```
   gronlund_driver: kmeans_hirschberg_larmore.cpp:67: virtual std::unique_ptr<kmeans_result> kmeans_hirschberg_larmore::compute(size_t): Assertion `intercept_guess > 0' failed.
   ```

   The last log line before the abort on N = 1000 was `call basic lambda=9.06305e-19` (the method prints `call basic lambda=...` to stdout on every λ trial).

2. **`wilber-binary` on `ties_n1000` k = 64 and `ties_n10000` k = 64** — still running at the 60 s cap, then SIGKILL (`timeout --signal=KILL`, rc = 137). GNU time: "Command terminated by signal 9", elapsed 1:00.00, peak RSS 1864 kB (N = 1000) and 1712 kB (N = 10,000). Stdout was only `setting search strategy to: 0`. No result was written.

3. **`wilber-interp` on `ties_n10000` k = 64 only** (N = 1000 succeeded) — SIGABRT, rc = 134, signal 6, elapsed 0.14 s. Stderr ended with:

   ```
   [Warning: K Found Outside search range - Empty Lambda Interval or numerical issues - Fix It with noise]
   stats [k_found, k-range searched] 404 - (20 ,  404 ) - lambda 1.1505555970355348749e-11
   gronlund_driver: kmeans_wilber.cpp:100: std::unique_ptr<kmeans_result> kmeans_wilber::compute_interpolation_search_with_noise(size_t, double): Assertion `succ' failed.
   ```

   Stdout had ten lines of `noise fail - try again - scale noise range by 1.5` after `setting search strategy to: 1`.

## Scale (N ∈ {1e5, 1e6}, k ∈ {64, 512})

One generated normal sample plus the audio prefixes. Every planned cell finished (64 cells, all `ok`). LEX is `lex_hybrid` with M = 4096. The five exact Grønlund methods agreed with each other on SSE to the printed digits in every cell (`hirc` / `wilber-*` sometimes differ in the last reported unit in the 17th significant digit, e.g. 0.9209612437794306 vs 0.9209612437794307).

`elapsed s` = GNU time wall clock. `driver s` = in-process compute time. Peak RSS = GNU time max RSS in MiB.

| dataset | N | k | method | status | elapsed s | driver s | peak RSS MiB | SSE |
|---|---:|---:|---|---|---:|---:|---:|---:|
| normal_n100000 | 100000 | 64 | LEX hybrid M=4096 | ok | 0.060 | — | 10.51 | 64.5105297266 |
| normal_n100000 | 100000 | 64 | linear compute (cost only) | ok | 0.560 | 0.544 | 12.84 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | linear report (Hirschberg reconstr.) | ok | 1.51 | 1.50 | 18.33 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | monotone compute (cost only) | ok | 0.480 | 0.462 | 8.27 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | monotone report (Hirschberg reconstr.) | ok | 1.33 | 1.32 | 13.34 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | hirc (Hirschberg–Larmore + λ) | ok | 0.310 | 0.287 | 8.66 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | wilber-binary | ok | 0.370 | 0.350 | 10.53 | 64.5095262195 |
| normal_n100000 | 100000 | 64 | wilber-interp | ok | 0.160 | 0.145 | 9.75 | 64.5095262195 |
| normal_n100000 | 100000 | 512 | LEX hybrid M=4096 | ok | 0.130 | — | 16.10 | 0.921451662737 |
| normal_n100000 | 100000 | 512 | linear compute (cost only) | ok | 4.92 | 4.90 | 12.99 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | linear report (Hirschberg reconstr.) | ok | 12.77 | 12.75 | 18.27 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | monotone compute (cost only) | ok | 4.75 | 4.73 | 8.29 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | monotone report (Hirschberg reconstr.) | ok | 11.03 | 11.00 | 13.34 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | hirc (Hirschberg–Larmore + λ) | ok | 0.440 | 0.423 | 8.61 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | wilber-binary | ok | 0.540 | 0.524 | 10.57 | 0.920961243779 |
| normal_n100000 | 100000 | 512 | wilber-interp | ok | 0.280 | 0.256 | 8.47 | 0.920961243779 |
| normal_n1000000 | 1000000 | 64 | LEX hybrid M=4096 | ok | 0.770 | — | 75.85 | 650.285461314 |
| normal_n1000000 | 1000000 | 64 | linear compute (cost only) | ok | 8.86 | 8.64 | 98.62 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | linear report (Hirschberg reconstr.) | ok | 17.06 | 16.83 | 152.03 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | monotone compute (cost only) | ok | 5.86 | 5.62 | 49.53 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | monotone report (Hirschberg reconstr.) | ok | 15.74 | 15.50 | 102.58 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | hirc (Hirschberg–Larmore + λ) | ok | 3.62 | 3.38 | 52.93 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | wilber-binary | ok | 3.44 | 3.20 | 75.33 | 650.282481626 |
| normal_n1000000 | 1000000 | 64 | wilber-interp | ok | 1.78 | 1.55 | 70.54 | 650.282481626 |
| normal_n1000000 | 1000000 | 512 | LEX hybrid M=4096 | ok | 0.800 | — | 75.90 | 10.151904481 |
| normal_n1000000 | 1000000 | 512 | linear compute (cost only) | ok | 54.77 | 54.55 | 98.72 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | linear report (Hirschberg reconstr.) | ok | 139.62 | 139.39 | 152.20 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | monotone compute (cost only) | ok | 43.87 | 43.64 | 49.47 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | monotone report (Hirschberg reconstr.) | ok | 129.55 | 129.31 | 102.39 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | hirc (Hirschberg–Larmore + λ) | ok | 2.44 | 2.16 | 52.93 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | wilber-binary | ok | 5.13 | 4.90 | 75.28 | 10.1469109777 |
| normal_n1000000 | 1000000 | 512 | wilber-interp | ok | 1.17 | 0.958 | 50.80 | 10.1469109777 |
| audio_prefix_n100000 | 100000 | 64 | LEX hybrid M=4096 | ok | 0.050 | — | 10.52 | 1.69157294911 |
| audio_prefix_n100000 | 100000 | 64 | linear compute (cost only) | ok | 0.560 | 0.547 | 12.92 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | linear report (Hirschberg reconstr.) | ok | 1.43 | 1.41 | 18.27 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | monotone compute (cost only) | ok | 0.490 | 0.472 | 8.27 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | monotone report (Hirschberg reconstr.) | ok | 1.29 | 1.27 | 13.29 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | hirc (Hirschberg–Larmore + λ) | ok | 0.300 | 0.289 | 8.59 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | wilber-binary | ok | 0.340 | 0.326 | 10.57 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 64 | wilber-interp | ok | 0.180 | 0.164 | 9.72 | 1.69157253487 |
| audio_prefix_n100000 | 100000 | 512 | LEX hybrid M=4096 | ok | 0.090 | — | 16.09 | 0.0252665249339 |
| audio_prefix_n100000 | 100000 | 512 | linear compute (cost only) | ok | 4.54 | 4.52 | 12.99 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | linear report (Hirschberg reconstr.) | ok | 11.70 | 11.68 | 18.34 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | monotone compute (cost only) | ok | 3.76 | 3.74 | 8.33 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | monotone report (Hirschberg reconstr.) | ok | 11.07 | 11.05 | 13.34 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | hirc (Hirschberg–Larmore + λ) | ok | 0.400 | 0.390 | 8.59 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | wilber-binary | ok | 0.700 | 0.675 | 10.53 | 0.0252558419064 |
| audio_prefix_n100000 | 100000 | 512 | wilber-interp | ok | 0.230 | 0.216 | 8.58 | 0.0252558419064 |
| audio_prefix_n1000000 | 1000000 | 64 | LEX hybrid M=4096 | ok | 0.640 | — | 75.47 | 22.2450381384 |
| audio_prefix_n1000000 | 1000000 | 64 | linear compute (cost only) | ok | 6.35 | 6.14 | 98.55 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | linear report (Hirschberg reconstr.) | ok | 15.87 | 15.65 | 152.17 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | monotone compute (cost only) | ok | 5.77 | 5.60 | 49.53 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | monotone report (Hirschberg reconstr.) | ok | 15.56 | 15.37 | 102.65 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | hirc (Hirschberg–Larmore + λ) | ok | 3.36 | 3.17 | 53.03 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | wilber-binary | ok | 3.63 | 3.44 | 75.39 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 64 | wilber-interp | ok | 1.50 | 1.32 | 62.65 | 22.2447803343 |
| audio_prefix_n1000000 | 1000000 | 512 | LEX hybrid M=4096 | ok | 0.610 | — | 75.47 | 0.348386029774 |
| audio_prefix_n1000000 | 1000000 | 512 | linear compute (cost only) | ok | 52.75 | 52.54 | 98.72 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | linear report (Hirschberg reconstr.) | ok | 124.53 | 124.36 | 152.19 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | monotone compute (cost only) | ok | 41.97 | 41.79 | 49.47 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | monotone report (Hirschberg reconstr.) | ok | 118.59 | 118.40 | 102.49 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | hirc (Hirschberg–Larmore + λ) | ok | 4.83 | 4.63 | 52.97 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | wilber-binary | ok | 5.49 | 5.31 | 75.28 | 0.34824943667 |
| audio_prefix_n1000000 | 1000000 | 512 | wilber-interp | ok | 2.46 | 2.25 | 50.97 | 0.34824943667 |

The full 8,168,256-row audio file was not run.

## Observations (numbers only)

- The authors' tree at `585386d` built with stock `make timing` on g++ 14.2.0, C++11, `-O2`. No source patch was needed.
- All of the table-free / O(n)-space variants listed above exist in this repo and ran to completion on the continuous (normal and audio-prefix) scale inputs, including N = 1e6 and k = 512, and they agreed with each other on SSE.
- On correctness cells with exact SSE > 1e-20, max |relative SSE vs LEX exact| = 4.58 × 10⁻¹¹.
- Failures were confined to the tie-heavy sets with k = 64 > 20 distinct values: `hirc` aborted on both N; `wilber-binary` was killed at 60 s on both N (RSS < 2 MiB); `wilber-interp` aborted on N = 10,000 only. The assert messages are quoted above.
- On those same k > D ties cells, the DP `report` methods that finished returned SSE = 0 and `n_centers = 23`.
- Peak RSS on N = 1e6 was about 50–152 MiB for the Grønlund methods and about 75–76 MiB for LEX hybrid M = 4096 (GNU time, includes loading the text file).
- At N = 1e6, k = 512, GNU-time elapsed was 0.80 s (normal) / 0.61 s (audio prefix) for LEX hybrid M = 4096, 0.96–2.25 s driver time for `wilber-interp`, 2.16–4.63 s for `hirc`, 4.90–5.31 s for `wilber-binary`, 43.6–54.6 s for `linear`/`monotone` cost-only, and 118–139 s for their Hirschberg `report` reconstructions.

## Files

- `kmeans1d/` — clone at the commit above
- `gronlund_driver`, `driver.cpp` — bench driver
- `lex_build/lex_cli` — paper LEX CLI
- `data/` — inputs
- `logs/` — per-cell stdout, stderr, GNU time, and driver result files
- `results/correctness.csv`, `results/scale.csv`
- `env.txt`, `build_timing.log`, `run_bench.py`, `gen_data.py`
