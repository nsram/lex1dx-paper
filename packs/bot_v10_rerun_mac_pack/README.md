# v10 rerun pack: scaling grid + controls, with raw-data Refine and k-preserving sig5

One pack that reruns the v8 270-cell scaling grid and the v9 controls (A, B, C),
using the library-matching sig5 pipeline:

- **point-mass DP**: groups enter as (W, Wμ, Wμ²)
- **raw-data Refine**: after the weighted DP, map bounds to raw indices and
  evaluate moves on distinct raw values (the same atoms Tuned uses)
- **k-preserving grouping**: sig5 starts at 5 significant digits; if that
  yields D′ < k and the raw data has D ≥ k distinct values, grouping is
  retried at 6, 7, … digits (up to 15), then raw distinct-value atoms, so
  the result has min(k, D) bins. `sig` (digits used; 0 = raw / Tuned) is
  logged on every CELL line.

**sig5exact** is unchanged as a control: always 5-significant-digit grouping,
exact weighted wilber DP on all D′ groups, no k-fallback, no Collapse, no Refine.

Unzip this pack into its own folder (not inside the lex1d-paper repo). From
the pack folder:

```
bash run/all.sh
```

Or one part at a time:

```
bash run/scaling.sh
bash run/controls.sh
```

The scripts find `python3` (3.8+, with NumPy 1.17+) and a C++ compiler
themselves. Do not pass paths. The machine is kept awake with `caffeinate`
on macOS. Each part is **resumable**: the same command skips cells already
`ok` in its `out/*/runs.csv`.

Put the audio file in `data/` as before (`data/dharini-abhogi.txt.gz`). The
scripts generate the synthetic files from fixed seeds (hashes checked against
`data/synth/MANIFEST.sha256`; the large files are not in the zip). GPT-2
h17/h35 are used in the scaling grid if they are in `data/` or can be fetched.
See `data/README.md`. A synthetic or audio hash mismatch is a warning, not
an error.

## Expected runtime (M3 Max)

- `bash run/scaling.sh` — about **1 hour** (v8 270-cell grid: N = 2^16..2^23
  synthetic + audio + GPT-2 if present; k ∈ {64, 256, 512}; 1 warm-up + 3
  timed reps; Grønlund on the large N dominates)
- `bash run/controls.sh` — about **16 minutes** (A/B/C on 3 × 8M synthetic +
  audio)
- `bash run/all.sh` — about **1 hour 20 minutes**

A short check (smallest N, a couple of k, 1 rep):

```
bash run/all.sh --smoke
```

uses N = 2^16..2^18 for scaling (no audio/GPT-2) and N = 2^16 + a small fake
audio for controls, and should finish in a few minutes.

## Experiments

All timing is **post-sort**: every method gets the same std::sort-ed array;
sort time is the `sort_s` reference column. All SSE is scored on the raw
sorted points; REE is vs Grønlund SSE* (patched wilber, unit weights, as in
v8). M = max(4096, 32k).

**Scaling grid** (v8 cells): tuned, sig5, and Grønlund wilber-interp at
k ∈ {64, 256, 512} on synthetic uniform / normal / lognormal at
N = 2^16..2^23, plus audio and GPT-2 h17/h35 if present.

**Controls** (v9 cells), all on 3 × 8M synthetic + audio unless noted:

- **A. Baseline A/B** (k = 64, 256, 512). `wilber_pristine` = unmodified
  `vendor/kmeans1d` (upstream 585386d) in its own binary `pristine_cell`.
  `wilber` = the patched `bench/weighted_wilber` copy with unit weights.
- **B. sig5exact control** (k = 64, 256, 512). Always 5-significant-digit
  grouping, exact weighted wilber DP on all D′ groups as point masses;
  **no k-fallback, no Collapse, no Refine**. Default `sig5` (with fallback
  and raw Refine) is re-run alongside.
- **C. Collapse quantile** (normal 8M and audio, k = 512). tuned and sig5 at
  q = 0.5, 0.75 (default), 0.9.

## Outputs

- `out/scaling/runs.csv`, `out/scaling/summary.md`, `out/scaling/run.log`
- `out/controls/runs.csv`, `out/controls/summary.md`, `out/controls/run.log`
- `out/env.txt` machine, compiler, Python, NumPy
- `out/scaling_mac_out.zip`, `out/controls_mac_out.zip` (per-part)
- `out/v10_mac_out.zip` (from `run/all.sh`; send this back)

Timeouts: Grønlund cells have a timeout (they can loop or crash on tie-heavy
or k>D data); a timeout or crash is recorded and the run continues.

lex1d freeze: commit `8a691ee005b657f985ffa299c57542b79570b08c`. Weighted
cost patch, tuned Collapse, raw-data Refine, and the k-preserving fallback
live only in `bench/`. Grønlund is vendored from `github.com/gronlund/kmeans1d`
at `585386d`.
