"""LEX: fast least-squares partition of 1-D data into k contiguous bins.

    import numpy as np, lex
    r = lex.lex_hybrid(np.random.default_rng(0).normal(size=10_000), 16)
    r.sizes, r.centers, r.sse, r.assign([0.0, 1.5])

``lex_hybrid`` (the default method) runs the exact DP on ``coarse_groups`` Collapse groups, then LEX
refinement. ``lex_plain`` is plain LEX (Collapse start, then refinement); the 2024 deck preset is
``lex_plain(x, k, precision=1e-4, multi_start=True)``. ``ckmeans_exact`` is the global optimum by
dynamic programming, as Ckmeans.1d.dp. All accept non-negative ``weights``. On ties among equal-SSE
partitions all return the largest cuts, as Ckmeans.1d.dp (R) does. A k above the number of distinct values (or of deck
groups) is clamped to that number, as Ckmeans.1d.dp does, with a ``UserWarning``; ``Result.k`` is
the number of bins used.
Errors (empty or non-finite input, k < 1, bad weights or options, or finite input too large for
double arithmetic: n * max |x - mean|, or with weights max(W, 1) * max |x - mean| and W, above
sqrt(DBL_MAX) / 4 = 3.35e153; the message gives the limit and the observed values) raise ``ValueError``.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import _core

__all__ = ["Result", "lex_hybrid", "lex_plain", "ckmeans_exact", "DEFAULT_COARSE_GROUPS"]
__version__ = "0.1.0"

DEFAULT_COARSE_GROUPS: int = _core.DEFAULT_COARSE_GROUPS

@dataclass(frozen=True)
class Result:
    """A fitted partition of the sorted data. Bin i holds sorted values bounds[i]:bounds[i+1]."""

    bounds: np.ndarray   # k + 1 cut indices into the sorted data (int64)
    sizes: np.ndarray    # values per bin (int64)
    centers: np.ndarray  # bin means (weighted means with weights)
    upper: np.ndarray    # largest value in each bin
    sse: float           # (weighted) sum of squared deviations from the bin means
    start: str           # collapse, coarse_exact, equal_frequency, equal_width or exact
    atoms: int           # distinct values (or groups when precision > 0)
    sweeps: int
    moves: int
    presorted: bool      # the input was already sorted, so the sort was skipped

    @property
    def k(self) -> int:
        """Bins used: min(requested k, atoms)."""
        return int(self.sizes.shape[0])

    def assign(self, x) -> np.ndarray:
        """Bin index for each value: the first bin whose upper value is >= x (last bin if none)."""
        idx = np.searchsorted(self.upper, np.asarray(x, dtype=np.float64), side="left")
        return np.minimum(idx, self.k - 1)


def _vec(a, name: str) -> np.ndarray:
    v = np.ascontiguousarray(a, dtype=np.float64)
    if v.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {v.shape}")
    return v


def _clamp_warn(r: Result, k, precision: float = 0.0) -> Result:
    """Warns (as Ckmeans.1d.dp does) when the core clamped k; called directly by the public functions."""
    if r.k < int(k):
        what = f"groups at precision {precision:g}" if precision > 0 else "distinct values"
        warnings.warn(f"k = {int(k)} exceeds the number of {what} ({r.atoms}); using k = {r.k}",
                      UserWarning, stacklevel=3)
    return r


def _lex(x, k, weights, precision, multi_start, tolerance, coarse_groups) -> Result:
    w = None if weights is None else _vec(weights, "weights")
    return Result(**_core.lex(_vec(x, "x"), int(k), w, float(precision), bool(multi_start), float(tolerance),
                              coarse_groups))


def lex_hybrid(x, k: int, *, weights=None, coarse_groups: int = DEFAULT_COARSE_GROUPS,
               precision: float = 0.0, multi_start: bool = False,
               tolerance: float = _core.DEFAULT_TOLERANCE) -> Result:
    """LEX with the coarse exact start: Collapse to ``coarse_groups`` groups, the exact DP on them,
    then refinement. ``coarse_groups`` <= k gives lex_plain's result."""
    if int(coarse_groups) < 0:
        raise ValueError(f"coarse_groups must be >= 0, got {coarse_groups}")
    return _clamp_warn(_lex(x, k, weights, precision, multi_start, tolerance, int(coarse_groups)), k, precision)


def lex_plain(x, k: int, *, weights=None, precision: float = 0.0, multi_start: bool = False,
              tolerance: float = _core.DEFAULT_TOLERANCE) -> Result:
    """Plain LEX: the Collapse start, then refinement."""
    return _clamp_warn(_lex(x, k, weights, precision, multi_start, tolerance, 0), k, precision)


def ckmeans_exact(x, k: int, *, weights=None) -> Result:
    """Globally optimal partition by dynamic programming (O(k d log d) after the sort)."""
    w = None if weights is None else _vec(weights, "weights")
    return _clamp_warn(Result(**_core.ckmeans_exact(_vec(x, "x"), int(k), w)), k)
