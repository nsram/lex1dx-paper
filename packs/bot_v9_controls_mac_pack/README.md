# v9 controls pack: baseline A/B, sig5exact control, Collapse quantile

Three control experiments for the lex1dx paper. Built from the v8 scaling pack
(same fixed sig5 driver, same data handling, same post-sort timing).

From the pack folder:

```
bash run/controls.sh
```

Put the audio file in `data/` as before (`data/dharini-abhogi.txt.gz`). The
script finds `python3` (3.8+, NumPy 1.17+) and a C++ compiler itself and
generates the 8M synthetic files (uniform / normal / lognormal, n = 8388608).
A synthetic or audio hash mismatch is a warning, not an error. No GPT-2.

Expected runtime: **about 20–25 minutes** on the M3 Max. Grønlund dominates:
v8 measured ~8–12 s per wilber run at 8M, 4 runs per cell (1 warm-up + 3 timed),
times 2 builds x 12 cells. Everything else is seconds. Resumable: re-running
the same command skips cells already `ok` in `out/controls/runs.csv`.

Quick check (N = 2^16 + a small fake audio, 1 rep, ~20 s):

```
bash run/controls.sh --smoke
```

## Experiments (all on 3 x 8M synthetic + audio unless noted)

All timing is **post-sort**: every method gets the same std::sort-ed array;
sort time is the `sort_s` reference column. All SSE is computed on the raw
sorted points; REE is vs Grønlund SSE* (patched wilber, unit weights, as in v8).

- **A. Baseline A/B** (k = 64, 256, 512). `wilber_pristine` = unmodified
  `vendor/kmeans1d` (upstream 585386d) built into its own binary
  `pristine_cell`. `wilber` = the patched `bench/weighted_wilber` copy with
  unit weights, which v7/v8 used as the baseline. Same driver source, same
  compiler flags. Reports both post_s, the ratio, and whether SSE is identical.
- **B. sig5exact control** (k = 64, 256, 512). 5-significant-digit grouping,
  then exact weighted wilber DP on all D′ groups as point masses
  (w, w·mu, w·mu²); **no Collapse, no Refine**. Cuts are mapped back to raw
  points and scored there. Default `sig5` is re-run alongside, and v8's sig5
  post_s (from `ref/v8_runs.csv`) is shown for comparison.
- **C. Collapse quantile** (normal 8M and audio, k = 512). tuned and sig5 at
  q = 0.5, 0.75 (default), 0.9. Reports post_s, REE, Collapse passes, G.

## Outputs

- `out/controls/runs.csv` one row per cell (dataset, method, k, q, experiment,
  sort_s, post_s, stage times, D, D′, G, Collapse passes/merges, Refine
  sweeps/moves, sse, exact_sse, ree_pct, peak_rss_mb, status, note)
- `out/controls/summary.md` tables A, B, C and a bit-identity check vs v8
- `out/controls/run.log` raw driver output and the command for every cell
- `out/env.txt` machine, compiler, Python, NumPy
- `out/controls_mac_out.zip` all of the above. Send this back.

Timeouts: none expected (limit 300 s per Grønlund cell at 8M, v8 took ≤ 48 s).
A timeout or crash is logged as `status=timeout`/`crash` with its limit, not
truncated; the run continues.

## Code changes vs v8 (pack only)

- `bench/tuned_collapse.hpp`: quantile is a parameter, default 0.75
  (bit-identical default behaviour).
- `bench/scale_cell.cpp`: method `sig5exact`; optional `--q=Q`; `q` in the
  CELL line; `-DLEX_PRISTINE` build that only runs wilber against
  `vendor/kmeans1d` and reports `wilber_pristine`.
- `CMakeLists.txt`: adds the `pristine_cell` target.
- `run/controls.sh`, `run/controls_run.py` (reuses `run/scaling_run.py`).

lex1d freeze: commit `8a691ee005b657f985ffa299c57542b79570b08c` (checked by hash).
