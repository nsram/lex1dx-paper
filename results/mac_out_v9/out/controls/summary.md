# v9 controls

Post-sort timing only (std::sort excluded, shown as sort_s). Times are medians of the timed reps. All SSE is on the raw sorted points; REE is vs Grønlund SSE* (patched wilber, unit weights).

Cells: 64 rows, 0 not ok.

## A. Baseline A/B: pristine vs patched Grønlund wilber-interp

| dataset | k | sort_s | pristine post_s | patched post_s | patched/pristine | SSE identical |
|---|---:|---:|---:|---:|---:|:---:|
| uniform_n8388608 | 64 | 0.529 | 7.692 | 7.838 | 1.019 | NO 170.65890856765512 vs 170.65890856764022 |
| uniform_n8388608 | 256 | 0.529 | 11.016 | 11.123 | 1.010 | yes |
| uniform_n8388608 | 512 | 0.529 | 11.955 | 12.057 | 1.009 | yes |
| normal_n8388608 | 64 | 0.569 | 8.439 | 8.571 | 1.016 | yes |
| normal_n8388608 | 256 | 0.569 | 5.703 | 5.880 | 1.031 | yes |
| normal_n8388608 | 512 | 0.569 | 9.175 | 9.227 | 1.006 | yes |
| lognormal_n8388608 | 64 | 0.565 | 8.961 | 9.226 | 1.030 | yes |
| lognormal_n8388608 | 256 | 0.565 | 10.448 | 10.638 | 1.018 | yes |
| lognormal_n8388608 | 512 | 0.565 | 8.860 | 8.910 | 1.006 | yes |
| audio | 64 | 0.500 | 9.306 | 9.477 | 1.018 | yes |
| audio | 256 | 0.500 | 10.867 | 11.046 | 1.016 | yes |
| audio | 512 | 0.500 | 11.600 | 11.660 | 1.005 | yes |

## B. sig5exact control (5-sig rounding + exact weighted DP on all D′ groups; no Collapse, no Refine)

| dataset | k | D′ | sig5exact post_s | sig5exact REE % | sig5exact RSS MB | sig5 post_s | sig5 REE % | sig5 RSS MB | v8 sig5 post_s | wilber post_s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uniform_n8388608 | 64 | 239077 | 0.308 | 0.00001 | 113.4 | 0.054 | 0.00034 | 98.6 | 0.056 | 7.838 |
| uniform_n8388608 | 256 | 239077 | 0.362 | 0.00020 | 109.5 | 0.072 | 0.00453 | 98.6 | 0.073 | 11.123 |
| uniform_n8388608 | 512 | 239077 | 0.369 | 0.00080 | 109.5 | 0.094 | 0.00358 | 98.7 | 0.093 | 12.057 |
| normal_n8388608 | 64 | 462310 | 0.556 | 0.00001 | 242.3 | 0.081 | 0.00014 | 129.2 | 0.083 | 8.571 |
| normal_n8388608 | 256 | 462310 | 0.616 | 0.00024 | 221.9 | 0.108 | 0.00165 | 130.0 | 0.109 | 5.880 |
| normal_n8388608 | 512 | 462310 | 0.578 | 0.00089 | 214.5 | 0.125 | 0.00204 | 129.9 | 0.125 | 9.227 |
| lognormal_n8388608 | 64 | 244439 | 0.290 | 0.00001 | 114.5 | 0.061 | 0.00044 | 99.4 | 0.061 | 9.226 |
| lognormal_n8388608 | 256 | 244439 | 0.360 | 0.00015 | 108.8 | 0.077 | 0.00206 | 99.6 | 0.077 | 10.638 |
| lognormal_n8388608 | 512 | 244439 | 0.366 | 0.00070 | 108.8 | 0.093 | 0.00207 | 99.6 | 0.093 | 8.910 |
| audio | 64 | 552603 | 0.296 | 0.00001 | 164.6 | 0.091 | 0.00034 | 135.7 | 0.091 | 9.477 |
| audio | 256 | 552603 | 0.739 | 0.00008 | 163.2 | 0.118 | 0.00168 | 137.1 | 0.116 | 11.046 |
| audio | 512 | 552603 | 0.862 | 0.00032 | 355.5 | 0.135 | 0.00209 | 137.1 | 0.135 | 11.660 |

## C. Collapse quantile sensitivity (k = 512; 0.75 is the default)

| dataset | method | q | post_s | collapse_s | REE % | collapse passes | G |
|---|---|---:|---:|---:|---:|---:|---:|
| normal_n8388608 | tuned | 0.5 | 0.499 | 0.367 | 0.00167 | 17 | 16384 |
| normal_n8388608 | tuned | 0.75 | 0.374 | 0.240 | 0.00212 | 12 | 16384 |
| normal_n8388608 | tuned | 0.9 | 0.341 | 0.209 | 0.00111 | 10 | 16384 |
| normal_n8388608 | sig5 | 0.5 | 0.126 | 0.022 | 0.00320 | 10 | 16384 |
| normal_n8388608 | sig5 | 0.75 | 0.125 | 0.014 | 0.00204 | 7 | 16384 |
| normal_n8388608 | sig5 | 0.9 | 0.106 | 0.010 | 0.00321 | 6 | 16384 |
| audio | tuned | 0.5 | 0.442 | 0.304 | 0.00153 | 16 | 16384 |
| audio | tuned | 0.75 | 0.353 | 0.218 | 0.00216 | 12 | 16384 |
| audio | tuned | 0.9 | 0.313 | 0.183 | 0.00260 | 10 | 16384 |
| audio | sig5 | 0.5 | 0.142 | 0.024 | 0.00357 | 10 | 16384 |
| audio | sig5 | 0.75 | 0.135 | 0.017 | 0.00209 | 7 | 16384 |
| audio | sig5 | 0.9 | 0.163 | 0.012 | 0.00196 | 6 | 16384 |

## Reproducibility vs v8 (same dataset, method, k, default q)

28/28 cells have final SSE bit-identical to v8 (max relative diff 0.00e+00). Expected all identical on the same Mac; on another platform last-digit diffs come from long double width in the raw-SSE sum.

Wall time: 16.2 min. See runs.csv (one row per cell), run.log (raw driver output), ../env.txt.
