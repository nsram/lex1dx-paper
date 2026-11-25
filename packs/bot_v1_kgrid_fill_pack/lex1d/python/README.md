# lex1d (Python bindings for LEX)

Python bindings (nanobind) for the LEX C++ core: fast least-squares partitioning of 1-D data into
k contiguous bins, plus the exact dynamic-programming optimum. See the repository README for the
algorithm, accuracy and speed numbers.

```python
import numpy as np, lex
x = np.random.default_rng(0).lognormal(size=100_000)
r = lex.lex_hybrid(x, 16)              # the default: exact DP on 4096 coarse groups, then LEX
p = lex.lex_plain(x, 16)               # plain LEX (Collapse start)
e = lex.ckmeans_exact(x, 16)           # global optimum (DP), as Ckmeans.1d.dp
w = lex.lex_hybrid(x, 16, weights=np.ones_like(x))
labels = r.assign(x)
```

Build from a full checkout of https://github.com/nsram/lex1d (the extension compiles `../cpp`):
`pip install ./python`.

The 2024 deck preset is `lex.lex_plain(x, k, precision=1e-4, multi_start=True)`.
`lex_hybrid` also takes `coarse_groups`; both LEX methods take `precision`, `multi_start` and
`tolerance`.

A k above the number of distinct values (or deck groups) is clamped to it with a `UserWarning`, as
Ckmeans.1d.dp does; `Result.k` is the number of bins used.
