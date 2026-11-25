"""Per-row 4-bit quantization of a synthetic embedding-style matrix. Generates its own data.

    python examples/embeddings.py

A 2,000 x 64 matrix whose rows come from three contrasting families (uniform, normal, and
lognormal-magnitude with random signs, i.e. heavy-tailed). Each row is quantized on its own to
16 levels (4 bits) with lex_hybrid, lex_plain, ckmeans_exact and a uniform min-max grid; the table reports
the mean squared error per family and in total, and the time.
"""
import time

import numpy as np

import lex

ROWS, DIM, K = 2_000, 64, 16


def matrix(rng):
    fam = rng.integers(0, 3, ROWS)
    m = np.empty((ROWS, DIM))
    m[fam == 0] = rng.uniform(-1, 1, ((fam == 0).sum(), DIM))
    m[fam == 1] = rng.normal(0, 0.5, ((fam == 1).sum(), DIM))
    n2 = (fam == 2).sum()
    m[fam == 2] = rng.choice([-1.0, 1.0], (n2, DIM)) * rng.lognormal(-1, 1, (n2, DIM))
    return m, fam


def uniform_rows(m, k):
    lo, hi = m.min(1, keepdims=True), m.max(1, keepdims=True)
    step = (hi - lo) / k
    q = lo + (np.minimum(np.floor((m - lo) / step), k - 1) + 0.5) * step
    return q


def quantize_rows(m, fit):
    q = np.empty_like(m)
    for i, row in enumerate(m):
        r = fit(row)
        q[i] = r.centers[r.assign(row)]
    return q


if __name__ == "__main__":
    m, fam = matrix(np.random.default_rng(64))
    methods = {
        "lex_hybrid": lambda row: lex.lex_hybrid(row, K),
        "lex_plain": lambda row: lex.lex_plain(row, K),
        "ckmeans_exact": lambda row: lex.ckmeans_exact(row, K),
    }
    out = {}
    for name, f in methods.items():
        t = time.perf_counter()
        out[name] = (quantize_rows(m, f), time.perf_counter() - t)
    t = time.perf_counter()
    out["uniform min-max"] = (uniform_rows(m, K), time.perf_counter() - t)
    names = ("uniform", "normal", "heavy-tailed")
    print(f"{ROWS} x {DIM} matrix, per-row k = {K} (4 bits); MSE per row family")
    print(f"{'method':<16}" + "".join(f"{n:>14}" for n in names) + f"{'all':>12}{'vs exact':>10}{'time':>9}")
    exact_all = np.mean((m - out["ckmeans_exact"][0]) ** 2)
    for name, (q, t) in out.items():
        err = (m - q) ** 2
        cols = "".join(f"{err[fam == j].mean():>14.3e}" for j in range(3))
        print(f"{name:<16}{cols}{err.mean():>12.3e}{err.mean() / exact_all:>9.3f}x{t:>8.2f}s")
