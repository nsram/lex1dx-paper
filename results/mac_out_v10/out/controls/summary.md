# v9 controls

Post-sort timing only (std::sort excluded, shown as sort_s). Times are medians of the timed reps. All SSE is on the raw sorted points; REE is vs Grønlund SSE* (patched wilber, unit weights).

Cells: 64 rows, 0 not ok.

## A. Baseline A/B: pristine vs patched Grønlund wilber-interp

| dataset | k | sort_s | pristine post_s | patched post_s | patched/pristine | SSE identical |
|---|---:|---:|---:|---:|---:|:---:|
| uniform_n8388608 | 64 | 0.566 | 7.701 | 7.847 | 1.019 | NO 170.65890856765512 vs 170.65890856764022 |
| uniform_n8388608 | 256 | 0.566 | 10.981 | 11.019 | 1.003 | yes |
| uniform_n8388608 | 512 | 0.566 | 11.953 | 12.000 | 1.004 | yes |
| normal_n8388608 | 64 | 0.572 | 8.548 | 8.614 | 1.008 | yes |
| normal_n8388608 | 256 | 0.572 | 5.755 | 5.704 | 0.991 | yes |
| normal_n8388608 | 512 | 0.572 | 9.182 | 9.173 | 0.999 | yes |
| lognormal_n8388608 | 64 | 0.565 | 8.941 | 9.184 | 1.027 | yes |
| lognormal_n8388608 | 256 | 0.565 | 10.227 | 10.586 | 1.035 | yes |
| lognormal_n8388608 | 512 | 0.565 | 8.831 | 8.875 | 1.005 | yes |
| audio | 64 | 0.501 | 9.350 | 9.450 | 1.011 | yes |
| audio | 256 | 0.501 | 10.880 | 10.801 | 0.993 | yes |
| audio | 512 | 0.501 | 11.593 | 11.605 | 1.001 | yes |

## B. sig5exact control (5-sig rounding + exact weighted DP on all D′ groups; no Collapse, no Refine)

| dataset | k | D′ | sig5exact post_s | sig5exact REE % | sig5exact RSS MB | sig5 post_s | sig5 REE % | sig5 RSS MB | v8 sig5 post_s | wilber post_s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uniform_n8388608 | 64 | 239077 | 0.310 | 0.00001 | 113.4 | 0.117 | 0.00024 | 419.0 | 0.056 | 7.847 |
| uniform_n8388608 | 256 | 239077 | 0.363 | 0.00020 | 109.5 | 0.145 | 0.00371 | 419.9 | 0.073 | 11.019 |
| uniform_n8388608 | 512 | 239077 | 0.367 | 0.00080 | 109.5 | 0.164 | 0.00130 | 445.5 | 0.093 | 12.000 |
| normal_n8388608 | 64 | 462310 | 0.556 | 0.00001 | 159.8 | 0.142 | 0.00010 | 449.6 | 0.083 | 8.614 |
| normal_n8388608 | 256 | 462310 | 0.613 | 0.00024 | 150.5 | 0.173 | 0.00106 | 451.0 | 0.109 | 5.704 |
| normal_n8388608 | 512 | 462310 | 0.573 | 0.00089 | 153.4 | 0.199 | 0.00132 | 552.1 | 0.125 | 9.173 |
| lognormal_n8388608 | 64 | 244439 | 0.290 | 0.00001 | 114.5 | 0.145 | 0.00045 | 421.4 | 0.061 | 9.184 |
| lognormal_n8388608 | 256 | 244439 | 0.360 | 0.00015 | 108.8 | 0.184 | 0.00156 | 423.7 | 0.077 | 10.586 |
| lognormal_n8388608 | 512 | 244439 | 0.363 | 0.00070 | 145.0 | 0.186 | 0.00150 | 422.0 | 0.093 | 8.875 |
| audio | 64 | 552603 | 0.296 | 0.00001 | 164.6 | 0.146 | 0.00034 | 488.1 | 0.091 | 9.450 |
| audio | 256 | 552603 | 0.739 | 0.00008 | 155.1 | 0.184 | 0.00184 | 450.1 | 0.116 | 10.801 |
| audio | 512 | 552603 | 0.851 | 0.00032 | 154.9 | 0.198 | 0.00191 | 524.5 | 0.135 | 11.605 |

## C. Collapse quantile sensitivity (k = 512; 0.75 is the default)

| dataset | method | q | post_s | collapse_s | REE % | collapse passes | G |
|---|---|---:|---:|---:|---:|---:|---:|
| normal_n8388608 | tuned | 0.5 | 0.501 | 0.369 | 0.00167 | 17 | 16384 |
| normal_n8388608 | tuned | 0.75 | 0.374 | 0.239 | 0.00212 | 12 | 16384 |
| normal_n8388608 | tuned | 0.9 | 0.339 | 0.207 | 0.00111 | 10 | 16384 |
| normal_n8388608 | sig5 | 0.5 | 0.196 | 0.022 | 0.00225 | 10 | 16384 |
| normal_n8388608 | sig5 | 0.75 | 0.197 | 0.014 | 0.00132 | 7 | 16384 |
| normal_n8388608 | sig5 | 0.9 | 0.196 | 0.010 | 0.00217 | 6 | 16384 |
| audio | tuned | 0.5 | 0.440 | 0.305 | 0.00153 | 16 | 16384 |
| audio | tuned | 0.75 | 0.350 | 0.217 | 0.00216 | 12 | 16384 |
| audio | tuned | 0.9 | 0.310 | 0.182 | 0.00260 | 10 | 16384 |
| audio | sig5 | 0.5 | 0.208 | 0.024 | 0.00315 | 10 | 16384 |
| audio | sig5 | 0.75 | 0.196 | 0.018 | 0.00191 | 7 | 16384 |
| audio | sig5 | 0.9 | 0.207 | 0.012 | 0.00182 | 6 | 16384 |

## Reproducibility vs v8 (same dataset, method, k, default q)

tuned and wilber should match v8 (same pipeline). sig5 is expected to differ: v10 Refines on raw distinct values and raises grouping precision when D' < k.
14/14 tuned/wilber cells have final SSE bit-identical to v8 (max relative diff 0.00e+00).

Wall time: 16.2 min. See runs.csv (one row per cell), run.log (raw driver output), ../env.txt.
