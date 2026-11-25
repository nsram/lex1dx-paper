# k-grid LEX fix: lex1dx sig5 (Refine on precision-group atoms)

One warm-up and 3 timed reps; times are medians, in seconds. **Comparison is post-sort only**: every method is given the same array, already sorted with `std::sort` (that sort is timed separately and logged as `sort_s`; it is not added to `post_s`). lex1dx **sig5** always groups to 5 significant digits and Refines on the precision-group atoms; bounds are mapped to raw indices only for SSE scoring (the reported configuration). M = max(4096, 32k): M stays 4096 for k ≤ 128 and grows with k above that, so these cells evaluate the default budget policy rather than isolating k. If an exact SSE* is supplied it is the two-pass long-double SSE of unweighted Grønlund wilber-interp on the sorted points, and REE is (SSE − SSE*) / SSE* × 100. This pack does not rerun the exact solver. Peak RSS is the cell process (`getrusage`).

This pack reruns only LEX sig5. Exact-solver SSE* and REE are left blank unless an exact value was passed in.

## Sort reference (`std::sort`, not in post-sort)

| dataset | n | D | sort s | peak RSS MiB |
|---|---:|---:|---:|---:|
| audio | 8168256 | 7295677 | 0.467 | 188.3 |
| lognormal_n1048576 | 1048576 | 1048576 | 0.057 | 25.2 |
| lognormal_n8388608 | 8388608 | 8388608 | 0.528 | 193.2 |
| normal_n1048576 | 1048576 | 1048576 | 0.057 | 25.2 |
| normal_n8388608 | 8388608 | 8388608 | 0.534 | 193.2 |
| uniform_n1048576 | 1048576 | 1048576 | 0.058 | 25.2 |
| uniform_n8388608 | 8388608 | 8388608 | 0.529 | 193.2 |

## Post-sort

| dataset | k | M | method | status | D | D' | group s | collapse s | DP s | refine s | post-sort s | vs wilber | SSE | REE % | passes | merges | refine moves | RSS MiB |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| lognormal_n1048576 | 16 | 4096 | sig5 | ok | 1048576 | 173142 | 0.013 | 0.006 | 0.003 | 0.000 | 0.025 |  | 80546.672609607966 |  | 6 | 169046 | 78 | 34.7 |
| lognormal_n1048576 | 32 | 4096 | sig5 | ok | 1048576 | 173142 | 0.013 | 0.006 | 0.001 | 0.000 | 0.023 |  | 19120.46730891378 |  | 6 | 169046 | 198 | 34.6 |
| lognormal_n1048576 | 128 | 4096 | sig5 | ok | 1048576 | 173142 | 0.013 | 0.006 | 0.008 | 0.000 | 0.031 |  | 1018.5198960189923 |  | 6 | 169046 | 858 | 61.0 |
| normal_n1048576 | 4 | 4096 | sig5 | ok | 1048576 | 281824 | 0.018 | 0.010 | 0.001 | 0.000 | 0.034 |  | 123209.53913612443 |  | 7 | 277728 | 103 | 54.0 |
| normal_n1048576 | 8 | 4096 | sig5 | ok | 1048576 | 281824 | 0.018 | 0.010 | 0.003 | 0.000 | 0.036 |  | 36261.230409299947 |  | 7 | 277728 | 238 | 54.1 |
| normal_n1048576 | 16 | 4096 | sig5 | ok | 1048576 | 281824 | 0.018 | 0.009 | 0.003 | 0.000 | 0.036 |  | 9977.0716197382553 |  | 7 | 277728 | 207 | 54.1 |
| normal_n1048576 | 32 | 4096 | sig5 | ok | 1048576 | 281824 | 0.018 | 0.010 | 0.003 | 0.000 | 0.036 |  | 2627.4048213562796 |  | 7 | 277728 | 572 | 54.1 |
| normal_n1048576 | 128 | 4096 | sig5 | ok | 1048576 | 281824 | 0.018 | 0.009 | 0.008 | 0.001 | 0.042 |  | 169.68920490504411 |  | 7 | 277728 | 1897 | 54.2 |
| uniform_n1048576 | 16 | 4096 | sig5 | ok | 1048576 | 158760 | 0.011 | 0.006 | 0.003 | 0.000 | 0.024 |  | 341.65782573454584 |  | 5 | 154664 | 201 | 32.3 |
| uniform_n1048576 | 32 | 4096 | sig5 | ok | 1048576 | 158760 | 0.011 | 0.006 | 0.005 | 0.000 | 0.026 |  | 85.389604424559963 |  | 5 | 154664 | 212 | 32.3 |
| uniform_n1048576 | 128 | 4096 | sig5 | ok | 1048576 | 158760 | 0.011 | 0.006 | 0.006 | 0.001 | 0.027 |  | 5.3231831803138476 |  | 5 | 154664 | 1219 | 32.4 |
| lognormal_n8388608 | 16 | 4096 | sig5 | ok | 8388608 | 244439 | 0.026 | 0.010 | 0.003 | 0.000 | 0.054 |  | 571018.05331781728 |  | 6 | 240343 | 178 | 99.3 |
| lognormal_n8388608 | 32 | 4096 | sig5 | ok | 8388608 | 244439 | 0.026 | 0.009 | 0.005 | 0.000 | 0.055 |  | 148299.17781918109 |  | 6 | 240343 | 279 | 99.4 |
| lognormal_n8388608 | 128 | 4096 | sig5 | ok | 8388608 | 244439 | 0.026 | 0.009 | 0.008 | 0.001 | 0.059 |  | 8982.9891059873207 |  | 6 | 240343 | 2313 | 99.4 |
| normal_n8388608 | 16 | 4096 | sig5 | ok | 8388608 | 462310 | 0.039 | 0.014 | 0.003 | 0.001 | 0.075 |  | 79694.922826514914 |  | 8 | 458214 | 1733 | 129.2 |
| normal_n8388608 | 32 | 4096 | sig5 | ok | 8388608 | 462310 | 0.039 | 0.014 | 0.004 | 0.001 | 0.077 |  | 20999.668109192229 |  | 8 | 458214 | 1768 | 129.3 |
| normal_n8388608 | 128 | 4096 | sig5 | ok | 8388608 | 462310 | 0.039 | 0.014 | 0.007 | 0.007 | 0.087 |  | 1367.6339276585725 |  | 8 | 458214 | 7451 | 183.4 |
| uniform_n8388608 | 16 | 4096 | sig5 | ok | 8388608 | 239077 | 0.025 | 0.008 | 0.003 | 0.000 | 0.051 |  | 2730.5683529021026 |  | 6 | 234981 | 170 | 98.6 |
| uniform_n8388608 | 32 | 4096 | sig5 | ok | 8388608 | 239077 | 0.025 | 0.008 | 0.005 | 0.000 | 0.052 |  | 682.52720568870689 |  | 6 | 234981 | 466 | 98.6 |
| uniform_n8388608 | 128 | 4096 | sig5 | ok | 8388608 | 239077 | 0.025 | 0.008 | 0.007 | 0.007 | 0.062 |  | 42.652820512237646 |  | 6 | 234981 | 2993 | 99.0 |
| audio | 16 | 4096 | sig5 | ok | 7295677 | 552603 | 0.045 | 0.017 | 0.003 | 0.000 | 0.084 |  | 2510.1000103370216 |  | 8 | 548507 | 824 | 135.7 |
| audio | 32 | 4096 | sig5 | ok | 7295677 | 552603 | 0.045 | 0.017 | 0.004 | 0.000 | 0.086 |  | 642.29272505037761 |  | 8 | 548507 | 879 | 222.6 |
| audio | 128 | 4096 | sig5 | ok | 7295677 | 552603 | 0.045 | 0.017 | 0.008 | 0.005 | 0.092 |  | 41.67378888732781 |  | 8 | 548507 | 4627 | 135.9 |

See `out/lex/runs.csv` and `out/lex/run.log`. Re-run the same command to resume; finished (status=ok) cells are skipped.
