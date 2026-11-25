# k-grid LEX fix: sig5 on the reported Refine path

Reruns only the LEX sig5 cells from the k-grid fill, using the reported
code path: group to 5 significant digits, Refine on the precision-group
atoms, and map bounds to raw indices only for SSE scoring. This is the
unmodified `scale_cell.cpp` (and `tuned_collapse.hpp`) from the harness
that produced the paper's LEX numbers. It does not Refine on raw
full-precision atoms, and it does not rerun wilber-interp or
Ckmeans.1d.dp.

Unzip this pack into its own folder. From the pack folder:

```
bash run/kgrid_lex_fix.sh
```

The script finds `python3` (3.8+, with NumPy 1.17+) and a C++ compiler
itself. Do not pass paths. The machine is kept awake with `caffeinate` on
macOS. The run is **resumable**: the same command skips cells already
`ok` or `skip_k_ge_D` in `out/lex/runs.csv`.

Put the audio file in `data/` as `data/dharini-abhogi.txt.gz`. The script
generates the synthetic files from fixed seeds (hashes checked against
`data/synth/MANIFEST.sha256`; the large files are not in the zip). See
`data/README.md`. A synthetic or audio hash mismatch is a warning, not an
error. If the audio file is missing, those cells are skipped and the
synthetic cells still run.

A short check (tiny N, 1 timed rep):

```
bash run/kgrid_lex_fix.sh --smoke
```

uses N = 4096 for all three families and should finish in a few minutes
the first time (compile) and in under a minute after that.

## Cells

| method | datasets | k |
|---|---|---|
| LEX (sig5) | uniform, normal, lognormal at N = 2^20 and 2^23; audio if present (N = 8,168,256) | 16, 32, 128 |
| LEX (sig5) | normal_n1048576 only | 4, 8 |

That is 23 LEX cells, plus one `std::sort` reference per dataset. Same
inputs, seeds, and timing protocol as the rest of the grid: post-sort,
one unrecorded warm-up and 3 timed repetitions, median (upper-middle
index). SSE is scored on the raw sorted points with the two-pass
long-double scorer. This pack does not compute exact-solver SSE*, so REE
is left blank.

M = max(4096, 32k). For every k in this pack, including 128, M stays at
4096. These results evaluate the default budget policy, not k in isolation.

## Expected runtime and memory (M3 Max, estimates)

Estimates are scaled from earlier sig5 runs on this code path, and are
labelled as estimates.

- Wall time: **about 5 minutes (estimate)**, including first-time
  generation of the 1M and 8M files. The 23 fits themselves are a small
  fraction of that.
- Peak memory: **about 130 MiB (estimate)**, from earlier 8M sig5 cells
  on this path (Refine works on the precision groups, not the raw N
  points).

## Outputs

- `out/env.txt` — machine, compiler, Python, NumPy, lex1d commit, source
  hashes for `lex.cpp` and the reported `scale_cell.cpp`, M policy,
  timed operations, seeds
- `out/lex/runs.csv`, `out/lex/summary.md`, `out/lex/run.log`
- `out/kgrid_lex_fix_mac_out.zip` — send this back

lex1d freeze: commit `8a691ee005b657f985ffa299c57542b79570b08c`. The
weighted cost patch and tuned Collapse live in `bench/`. Compiler flags
include `-ffp-contract=off`. Timed operation: the method call on an
already-sorted array.
