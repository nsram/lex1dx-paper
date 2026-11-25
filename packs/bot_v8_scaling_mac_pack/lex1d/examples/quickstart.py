"""LEX quickstart: three contrasting distributions side by side. Generates its own data.

    python examples/quickstart.py

For uniform, normal and lognormal samples (N = 200,000, 10 % of values rounded to 2 decimals so
there are ties), fits k = 16 bins with lex_hybrid (the default), lex_plain and ckmeans_exact, and prints:
the gap to the exact optimum, how much worse a uniform min-max grid is (MSE ratio), and time.
Then a weighted variant, and the cuts and centres for the normal sample.
"""
import time

import numpy as np

import lex

N, K = 200_000, 16


def sample(name, rng):
    x = {"uniform": lambda: rng.uniform(0, 1, N),
         "normal": lambda: rng.normal(0, 1, N),
         "lognormal": lambda: rng.lognormal(0, 1, N)}[name]()
    tie = rng.random(N) < 0.10
    x[tie] = np.round(x[tie], 2)  # some exact duplicates
    return x


def uniform_minmax_mse(x, k, w=None):
    """MSE of a k-level uniform grid over [min, max] (each value -> its cell's midpoint)."""
    lo, hi = x.min(), x.max()
    step = (hi - lo) / k
    q = lo + (np.minimum(np.floor((x - lo) / step), k - 1) + 0.5) * step
    return np.average((x - q) ** 2, weights=w)


def timed(f, *a, **kw):
    t = time.perf_counter()
    r = f(*a, **kw)
    return r, time.perf_counter() - t


def table(title, weighted):
    rng = np.random.default_rng(2026)
    print(title)
    print(f"{'dist':<10} {'N':>8} {'k':>3} {'gap hybrid':>12} {'gap plain':>10} {'uniform/LEX':>12} "
          f"{'t hybrid':>10} {'t plain':>8} {'t exact':>8}")
    for name in ("uniform", "normal", "lognormal"):
        x = sample(name, rng)
        w = rng.lognormal(0, 1, N) if weighted else None
        d, td = timed(lex.lex_hybrid, x, K, weights=w)
        p, tp = timed(lex.lex_plain, x, K, weights=w)
        e, te = timed(lex.ckmeans_exact, x, K, weights=w)
        wsum = N if w is None else w.sum()
        ratio = uniform_minmax_mse(x, K, w) / (d.sse / wsum)
        print(f"{name:<10} {N:>8} {K:>3} {100 * (d.sse / e.sse - 1):>11.2e}% {100 * (p.sse / e.sse - 1):>9.3f}% "
              f"{ratio:>11.2f}x {td * 1e3:>8.1f}ms {tp * 1e3:>6.1f}ms {te * 1e3:>6.1f}ms")
    print()


if __name__ == "__main__":
    lex.lex_hybrid(np.arange(100.0), 4)  # warm-up (first call loads the extension)
    table("unweighted", weighted=False)
    table("weighted (lognormal weights)", weighted=True)
    x = sample("normal", np.random.default_rng(7))
    r = lex.lex_hybrid(x, K)
    np.set_printoptions(precision=3, suppress=True, linewidth=100)
    print("normal sample, k = 16")
    print("cuts (upper value of each bin):", r.upper)
    print("centres:", r.centers)
    print("sizes:", r.sizes, " start:", r.start, " atoms:", r.atoms)
