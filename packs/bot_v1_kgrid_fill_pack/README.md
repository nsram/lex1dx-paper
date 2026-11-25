# k-grid fill: LEX, wilber-interp, and Ckmeans.1d.dp

Fills the missing k cells so the paper can report a k grid of
16, 32, 64, 128, 256, 512, plus a small-k check at 4 and 8.

This pack does not repeat cells that already exist at k = 64, 256, 512
(LEX and wilber-interp) or at k = 16, 64, 128, 256, 512 (Ckmeans.1d.dp).

Unzip this pack into its own folder. From the pack folder:

```
bash run/kgrid_fill.sh
```

The script finds `python3` (3.8+, with NumPy 1.17+), a C++ compiler, and
`Rscript` itself. It installs `Ckmeans.1d.dp` from CRAN into a pack-local
`rlib/` if it is missing. Do not pass paths. The machine is kept awake with
`caffeinate` on macOS. The run is **resumable**: the same command skips
cells already `ok` or `skip_k_ge_D` in `out/lex/runs.csv` and
`out/ckmeans/runs.csv`. A Ckmeans cell recorded as `skipped` (over a memory
or time cap) is retried if the caps now allow it.

Put the audio file in `data/` as `data/dharini-abhogi.txt.gz`. The script
generates the synthetic files from fixed seeds (hashes checked against
`data/synth/MANIFEST.sha256`; the large files are not in the zip). See
`data/README.md`. A synthetic or audio hash mismatch is a warning, not an
error. If the audio file is missing, those cells are skipped and the
synthetic cells still run.

A short check (tiny N, 1 timed rep):

```
bash run/kgrid_fill.sh --smoke
```

uses N = 4096 for all three families, runs LEX, wilber-interp, and
Ckmeans, and should finish in a few minutes the first time (compile, and
a CRAN install if needed) and in under a minute after that.

## Cells

LEX is the reported configuration: 5-significant-digit grouping, the
k-preserving digit fallback, and group-level Refine. The exact solver is
unweighted Grønlund wilber-interp (interpolation search). Ckmeans uses
`method = "linear"`.

| method | datasets | k | notes |
|---|---|---|---|
| LEX (sig5) and wilber-interp | uniform, normal, lognormal at N = 2^20 and 2^23; audio if present (N = 8,168,256) | 16, 32, 128 | post-sort, same inputs and timing as the existing scaling runs |
| Ckmeans.1d.dp | the same datasets | 32 | sorted and raw input |
| LEX, wilber-interp, and Ckmeans | normal_n1048576 only | 4, 8 | small-k check; Ckmeans on sorted and raw |

That is 46 LEX/wilber cells and 18 Ckmeans cells (64 method cells), plus
one `std::sort` reference per dataset. All timing is post-sort: every
method gets the same `std::sort`-ed array; sort time is the `sort_s`
column and is not added to `post_s`. One unrecorded warm-up and 3 timed
repetitions; the reported time is the median (upper-middle index). SSE is
scored on the raw sorted points with the two-pass long-double scorer.
REE is vs wilber-interp SSE*.

M = max(4096, 32k). For every k in this pack, including 128, M stays at
4096. These results evaluate the default budget policy, not k in isolation.

## Expected runtime and memory (M3 Max, estimates)

Estimates are scaled from earlier runs at nearby k and N, and are labelled
as estimates.

- Wall time: **about 30 minutes (estimate)**. LEX is a small fraction of
  that. Most of the time is wilber-interp at N ≈ 8e6 (about 10–15 minutes,
  estimate) and the Ckmeans k = 32 fits plus their wilber SSE* references
  (about 15 minutes, estimate). Add a few minutes the first time if
  `Ckmeans.1d.dp` has to be installed from CRAN.
- Peak memory: **about 4.2 GiB (estimate)** for Ckmeans k = 32 at
  N ≈ 8e6 (`RSS ≈ 0.22 GiB + 16 bytes × N × k`). LEX and wilber-interp
  peaked at about 1.1 GiB in earlier 8M cells. The default Ckmeans memory
  cap of 24 GiB is enough. You do not need to set `LEX_MEM_CAP_GIB`.

## Outputs

- `out/env.txt` — machine, compiler, Python, NumPy, R, package version and
  method, lex1d commit and source hash, M policy, timed operations, seeds
- `out/lex/runs.csv`, `out/lex/summary.md`, `out/lex/run.log`
- `out/ckmeans/runs.csv`, `out/ckmeans/summary.md`, `out/ckmeans/run.log`
- `out/kgrid_fill_mac_out.zip` — send this back

lex1d freeze: commit `8a691ee005b657f985ffa299c57542b79570b08c`. Weighted
cost patch, tuned Collapse, group-level Refine, and the k-preserving
fallback live in `bench/`. Grønlund is vendored from
`github.com/gronlund/kmeans1d` at `585386d`. Compiler flags include
`-ffp-contract=off`. Timed operation: the method call on an already-sorted
array (Ckmeans still sorts internally; that internal sort is not excluded).
