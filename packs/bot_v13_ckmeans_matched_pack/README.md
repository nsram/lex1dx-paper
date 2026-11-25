# v11 pack: Ckmeans.1d.dp on the lex1dx v10 inputs and protocol

A small pack that reruns `Ckmeans.1d.dp` (R) on cells matched to the v10
scaling study, so the archived Table 4 can be replaced with numbers that use
the same samples, the same timing protocol, and the same SSE scorer as the
rest of the paper.

Unzip this pack into its own folder. From the pack folder:

```
bash run/ckmeans_matched.sh
```

Or, sorted-input cells only (the v10-comparable post-sort times):

```
bash run/ckmeans_matched.sh --sorted-only
```

The scripts find `python3` (3.8+, with NumPy 1.17+), a C++ compiler, and
`Rscript` themselves. They install `Ckmeans.1d.dp` from CRAN into a pack-local
`rlib/` if it is missing. Do not pass paths. The machine is kept awake with
`caffeinate` on macOS. The run is **resumable**: the same command skips cells
already `ok`, `skipped`, or `skip_k_ge_D` in `out/ckmeans_matched/runs.csv`.

Put the audio file in `data/` as in v10 (`data/dharini-abhogi.txt.gz`). The
scripts generate the synthetic files from the v10 seeds (hashes checked
against `data/synth/MANIFEST.sha256`; the large files are not in the zip). See
`data/README.md`. A synthetic or audio hash mismatch is a warning, not an
error.

## Expected runtime (M3 Max)

- `bash run/ckmeans_matched.sh` — about **45 minutes** (3 families at
  N = 2^20 and 2^23, plus audio if present; k ∈ {16, 128, 512}; both
  pre-sorted and unsorted fits; 1 warm-up + 3 timed reps). Cells whose
  estimated peak RSS exceeds the memory cap are skipped (see below).
- `bash run/ckmeans_matched.sh --sorted-only` — about **25 minutes**.

A short check (tiny N, a couple of k, 1 rep):

```
bash run/ckmeans_matched.sh --smoke
```

uses N = 4096 for all three families, k ∈ {4, 16}, and should finish in a
few minutes the first time (CRAN install + compile) and in under a minute
after that.

## Cells and protocol

| dataset | N | k | inputs |
|---|---:|---|---|
| uniform, normal, lognormal | 1,048,576 (2^20) | 16, 128, 512 | sorted and raw |
| uniform, normal, lognormal | 8,388,608 (2^23) | 16, 128, 512 | sorted and raw |
| audio (if present) | 8,168,256 | 16, 128, 512 | sorted and raw |

That is 21 dataset×k cells, twice if both inputs are run, 42 in all. Under
the default 24 GiB cap, k = 512 at N ≈ 8 million is skipped (estimated
peak RSS about 64 GiB), leaving 18 dataset×k cells (36 fits with both inputs).

**Timing.** One unrecorded warm-up and 3 timed repetitions; the reported
time is the median (same upper-middle index as v10). Each cell is its own
R process. Peak RSS is that process (`wait4` / `getrusage`), in MiB.

**Sorted vs unsorted.** Every method in v10 was timed on an array already
sorted with `std::sort`; that sort is logged here as `sort_s` and is not
added to the fit time. `Ckmeans.1d.dp` always sorts internally (the docs
give O(n log n + kn), or O(kn) if `x` is already sorted). The internal sort
cannot be excluded, so the comparable **post-sort** number is a call on the
v10-sorted array (`input=sorted`). The unsorted number (`input=raw`) includes
the internal sort of the unsorted file, as in the archived Table 4.

**Method.** The pack uses `method = "linear"`, the package's fastest exact
solver (SMAWK-style, O(kn) after sort). It uses the most memory of the three
exact methods, still O(kn). The method and the installed package version are
written to `out/env.txt` and to every row. To override: `LEX_CKMEANS_METHOD=loglinear`.

**SSE and REE.** After each fit, bin edges are scored on the raw sorted
values with the same two-pass long-double `raw_sse` as v10 `scale_cell.cpp`.
Exact SSE* is unweighted Grønlund wilber-interp (interpolation search, the
patched v10 solver) scored with that same function, so REE is comparable to
the main-grid figures. Shipped `ref/v10_wilber_sse.csv` holds the Mac v10
SSE* for k ∈ {64, 256, 512} as a cross-check; this pack recomputes SSE* for
every k it runs, including 16 and 128, which v10 did not measure.

## Memory and time guards

`Ckmeans.1d.dp` with `method = "linear"` needed 1.2 to 8.3 GiB at N = 2^20
in the archived study, growing with k. At N = 2^23 that scales to well above
a typical laptop budget at k = 512. Before each cell the pack estimates:

    peak RSS ≈ 0.22 GiB + 16 bytes × N × k
    fit time ≈ 21.6 s × (N / 2^20) × (k / 512)

(the fit is to the archived 4.3.6 sorted runs at N = 2^20). If the estimate
exceeds a cap, the cell is not started and is recorded as
`status=skipped` with `note=skipped: reason`. Defaults:

```
LEX_MEM_CAP_GIB=24
LEX_TIME_CAP_S=1800
```

Raise the memory cap only if the machine has the RAM, for example
`LEX_MEM_CAP_GIB=80 bash run/ckmeans_matched.sh`. A process that still
overruns the estimated cap is timed out at `LEX_TIME_CAP_S` and recorded as
`status=timeout`.

## Outputs

- `out/ckmeans_matched/runs.csv`, `out/ckmeans_matched/summary.md`, `out/ckmeans_matched/run.log`
- `out/env.txt` machine, R version, package version and method, compiler, Python, NumPy, caps
- `out/ckmeans_matched_mac_out.zip` (send this back)
