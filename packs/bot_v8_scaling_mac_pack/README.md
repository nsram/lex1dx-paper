# Scaling pack: lex1dx vs Grønlund, N = 2^16..2^23

Post-sort comparison of:

- **lex1dx tuned** — full-precision atoms, tuned batch Collapse, weighted wilber-interp DP, Refine
- **lex1dx sig5** — always groups to 5 significant digits (no D′/M guard), then the same pipeline; D and D′ are logged in every cell so a rule can be derived later
- **Grønlund wilber-interp** — unweighted, on the raw sorted points (this SSE* is also the exact column)

M = max(4096, 32k). k ∈ {64, 256, 512}. Every method is given the **same pre-sorted** array. `std::sort` is timed separately and written as the `sort_s` reference column; it is not included in `post_s`. No radix sort.

Data: synthetic uniform / normal / lognormal at N = 2^16..2^23 (generated on the Mac from fixed seeds, hashes checked against `data/synth/MANIFEST.sha256`; the large files are not in the zip), plus audio and GPT-2 h17/h35 if they are in `data/` or can be fetched. See `data/README.md`.

Each cell records sort time, grouping time, D, D′, Collapse passes and merges per pass, Refine sweeps and moves, DP time, peak RSS, final SSE, and REE vs Grønlund SSE*. Grønlund cells have a timeout (they can loop or crash on tie-heavy or k>D data); a timeout or crash is recorded and the run continues.

lex1d freeze: commit `8a691ee005b657f985ffa299c57542b79570b08c`. Weighted cost patch and tuned Collapse live only in `bench/`. Grønlund is vendored from `github.com/gronlund/kmeans1d` at `585386d`.

Unzip this pack into its own folder (not inside the lex1d-paper repo). From the pack folder:

```
bash run/scaling.sh
```

The script finds `python3` (3.8+, with NumPy 1.17+) and a C++ compiler itself. Do not pass paths. The machine is kept awake with `caffeinate` on macOS. The run is **resumable**: the same command skips cells already in `out/scaling/runs.csv` with `status=ok`.

Expected runtime is about **25–40 minutes** on an M3 Max for the full set (24 synthetic files + audio + two GPT-2 tensors, 1 warm-up + 3 timed reps; wilber-interp on the large N dominates). Without audio and GPT-2, closer to **20–30 minutes**. A short check:

```
bash run/scaling.sh --smoke
```

uses N = 2^16..2^18 only, one warm-up and one timed rep, no audio/GPT-2, and should finish in a few minutes.

Results: `out/scaling/summary.md`, `out/scaling/runs.csv`, and `out/scaling_mac_out.zip`. Send that zip back.
