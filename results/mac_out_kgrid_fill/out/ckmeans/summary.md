# Ckmeans.1d.dp k-grid fill

One warm-up and 3 timed reps; times are medians, in seconds. Inputs, seeds, and files are those of the v10 scaling pack (`run/gen_synth.py`, SeedSequence(2048, spawn_key=(stream_id, N)); audio is the v10 Dharini Abhogi recording, N = 8,168,256). `Ckmeans.1d.dp` ran with `method = "linear"` (fastest exact method, O(kn) after sort, O(kn) memory). **Post-sort time (`input=sorted`)** is a call on an array already sorted with `std::sort` (timed separately as `sort_s`). `Ckmeans.1d.dp` still sorts internally; that internal sort is on already-sorted data and is **not excluded**. Unsorted fits (`input=raw`) include the internal sort of the unsorted array, matching the archived Table 4 column. SSE is scored on the raw sorted values with the same two-pass long-double `raw_sse` as v10 `scale_cell.cpp`. Exact SSE* is unweighted Grønlund wilber-interp (interpolation search, the v10 patched solver) scored with the same function. REE is (SSE − SSE*) / SSE* × 100. Peak RSS is the cell process (`wait4` / `getrusage`), in MiB. A cell is skipped, and recorded as `status=skipped`, if the estimated peak RSS exceeds 24 GiB (`LEX_MEM_CAP_GIB`) or the estimated wall time for the warm-up plus reps exceeds 1800 s (`LEX_TIME_CAP_S`). The estimator is a fit to the archived Ckmeans.1d.dp 4.3.6 runs at N = 2^20: RSS ≈ 0.22 GiB + 16 bytes × N × k, and time ≈ 21.6 s × (N / 2^20) × (k / 512) per fit.

## Sort reference (`std::sort`, not in post-sort)

| dataset | n | D | sort s | peak RSS MiB |
|---|---:|---:|---:|---:|
| audio | 8168256 | 7295677 | 0.501 | 188.2 |
| lognormal_n1048576 | 1048576 | 1048576 | 0.061 | 25.2 |
| lognormal_n8388608 | 8388608 | 8388608 | 0.567 | 193.2 |
| normal_n1048576 | 1048576 | 1048576 | 0.062 | 25.2 |
| normal_n8388608 | 8388608 | 8388608 | 0.577 | 193.3 |
| uniform_n1048576 | 1048576 | 1048576 | 0.062 | 25.2 |
| uniform_n8388608 | 8388608 | 8388608 | 0.569 | 193.2 |

## Ckmeans.1d.dp

| dataset | k | input | status | n | D | post-sort / fit s | SSE | SSE* | REE % | tot.withinss | RSS MiB | est. RSS GiB | version | note |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| lognormal_n1048576 | 32 | sorted | ok | 1048576 | 1048576 | 1.464 | 19120.463816152027 | 19120.463816152027 | 0 | 19120.46381615206 | 1057.3 | 0.72 | 4.3.6 |  |
| lognormal_n1048576 | 32 | raw | ok | 1048576 | 1048576 | 1.557 | 19120.463816152027 | 19120.463816152027 | 0 | 19120.46381615206 | 1102.8 | 0.72 | 4.3.6 |  |
| normal_n1048576 | 4 | sorted | ok | 1048576 | 1048576 | 0.125 | 123209.53873590003 | 123209.53873590003 | 0 | 123209.53873590007 | 510.8 | 0.28 | 4.3.6 |  |
| normal_n1048576 | 4 | raw | ok | 1048576 | 1048576 | 0.212 | 123209.53873590003 | 123209.53873590003 | 0 | 123209.53873590007 | 371.4 | 0.28 | 4.3.6 |  |
| normal_n1048576 | 8 | sorted | ok | 1048576 | 1048576 | 0.334 | 36261.230074745159 | 36261.230074745159 | 0 | 36261.230074744293 | 357.4 | 0.34 | 4.3.6 |  |
| normal_n1048576 | 8 | raw | ok | 1048576 | 1048576 | 0.388 | 36261.230074745159 | 36261.230074745159 | 0 | 36261.230074744293 | 493.4 | 0.34 | 4.3.6 |  |
| normal_n1048576 | 32 | sorted | ok | 1048576 | 1048576 | 1.419 | 2627.4041197454289 | 2627.4041197454289 | 0 | 2627.4041197454326 | 1397.7 | 0.72 | 4.3.6 |  |
| normal_n1048576 | 32 | raw | ok | 1048576 | 1048576 | 1.491 | 2627.4041197454289 | 2627.4041197454289 | 0 | 2627.4041197454326 | 1048.9 | 0.72 | 4.3.6 |  |
| uniform_n1048576 | 32 | sorted | ok | 1048576 | 1048576 | 1.361 | 85.389596173087924 | 85.389596173087924 | 0 | 85.389596173087469 | 964.3 | 0.72 | 4.3.6 |  |
| uniform_n1048576 | 32 | raw | ok | 1048576 | 1048576 | 1.442 | 85.389596173087924 | 85.389596173087924 | 0 | 85.389596173087469 | 915.8 | 0.72 | 4.3.6 |  |
| lognormal_n8388608 | 32 | sorted | ok | 8388608 | 8388608 | 11.808 | 148299.17639675189 | 148299.17639675189 | 0 | 148299.17639675547 | 5180.2 | 4.22 | 4.3.6 |  |
| lognormal_n8388608 | 32 | raw | ok | 8388608 | 8388608 | 12.842 | 148299.17639675189 | 148299.17639675189 | 0 | 148299.17639675547 | 5122.2 | 4.22 | 4.3.6 |  |
| normal_n8388608 | 32 | sorted | ok | 8388608 | 8388608 | 11.632 | 20999.64815871192 | 20999.64815871192 | 0 | 20999.648158713564 | 5223.3 | 4.22 | 4.3.6 |  |
| normal_n8388608 | 32 | raw | ok | 8388608 | 8388608 | 12.572 | 20999.64815871192 | 20999.64815871192 | 0 | 20999.648158713564 | 5116.1 | 4.22 | 4.3.6 |  |
| uniform_n8388608 | 32 | sorted | ok | 8388608 | 8388608 | 11.280 | 682.52505839127991 | 682.52505839127991 | 0 | 682.52505839127195 | 5208.0 | 4.22 | 4.3.6 |  |
| uniform_n8388608 | 32 | raw | ok | 8388608 | 8388608 | 12.346 | 682.52505839127991 | 682.52505839127991 | 0 | 682.52505839127195 | 5173.2 | 4.22 | 4.3.6 |  |
| audio | 32 | sorted | ok | 8168256 | 7295677 | 11.001 | 642.29266809745172 | 642.29266809745172 | 0 | 642.29266809742251 | 5049.3 | 4.12 | 4.3.6 |  |
| audio | 32 | raw | ok | 8168256 | 7295677 | 11.870 | 642.29266809745172 | 642.29266809745172 | 0 | 642.29266809742251 | 5004.3 | 4.12 | 4.3.6 |  |

See `out/ckmeans/runs.csv` and `out/ckmeans/run.log`. Re-run the same command to resume; finished (`ok`, `skip_k_ge_D`) cells are skipped; `skipped` cells are re-checked against the current caps.
