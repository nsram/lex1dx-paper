import warnings
import struct
import subprocess

import numpy as np
import pytest

import lex


def load(path):
    return np.array([float(s) for s in path.read_text().split()], dtype=np.float64)


def cases(fixtures_dir):
    for line in (fixtures_dir / "cases.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        name, k, w = line.split()
        yield name, int(k), (None if w == "-" else w)


# lex_cli preset -> the Python call it stands for
METHODS = {
    "hybrid": lex.lex_hybrid,
    "plain": lex.lex_plain,
    "deck": lambda x, k, **kw: lex.lex_plain(x, k, precision=1e-4, multi_start=True, **kw),
    "exact": lex.ckmeans_exact,
}


def bits(x):
    return format(struct.unpack("<Q", struct.pack("<d", x))[0], "016x")


def run_py(fixtures_dir, name, k, w, preset):
    x = load(fixtures_dir / "data" / f"{name}.txt")
    wv = None if w is None else load(fixtures_dir / "data" / f"{w}.txt")
    try:
        r = METHODS[preset](x, k, weights=wv)
    except ValueError as e:
        return {"error": str(e)}
    return {
        "bounds": ",".join(map(str, r.bounds)),
        "sizes": ",".join(map(str, r.sizes)),
        "sse_bits": bits(r.sse),
        "centers_bits": ",".join(bits(c) for c in r.centers),
        "start": r.start,
        "atoms": str(r.atoms),
        "sweeps": str(r.sweeps),
        "moves": str(r.moves),
    }


def parse(text):
    return dict(line.split("=", 1) for line in text.splitlines() if "=" in line)


def test_bit_identical_to_cpp_cli(fixtures_dir, lex_cli):
    n = 0
    for name, k, w in cases(fixtures_dir):
        for preset in METHODS:
            args = [str(lex_cli), str(fixtures_dir / "data" / f"{name}.txt"), str(k), "--preset", preset]
            if w:
                args += ["--weights", str(fixtures_dir / "data" / f"{w}.txt")]
            cpp = parse(subprocess.run(args, capture_output=True, text=True).stdout)
            cpp.pop("sse", None)
            assert run_py(fixtures_dir, name, k, w, preset) == cpp, (name, k, w, preset)
            n += 1
    assert n >= 352  # 88 cases x 4 presets, weighted and clamped cases included


def test_expected_file_is_current(fixtures_dir, lex_cli):
    """fixtures/expected_cpp.txt (read by the Scala cross-check) matches today's C++ CLI."""
    for line in (fixtures_dir / "expected_cpp.txt").read_text().splitlines():
        if line.startswith("#"):
            continue
        name, k, w, preset, *kv = line.split("|")
        args = [str(lex_cli), str(fixtures_dir / "data" / f"{name}.txt"), k, "--preset", preset]
        if w != "-":
            args += ["--weights", str(fixtures_dir / "data" / f"{w}.txt")]
        out = subprocess.run(args, capture_output=True, text=True).stdout.splitlines()
        # sse=%.17g is a readable copy (printf formatting may differ by C runtime); sse_bits is exact
        drop = lambda fields: [f for f in fields if not f.startswith("sse=")]
        assert drop(out) == drop(kv), (name, k, w, preset)


def test_agrees_with_ckmeans_exact(fixtures_dir):
    """lex_hybrid (coarse start) never beats the optimum and stays within 0.01 % of it here."""
    for name, k, w in cases(fixtures_dir):
        x = load(fixtures_dir / "data" / f"{name}.txt")
        wv = None if w is None else load(fixtures_dir / "data" / f"{w}.txt")
        e = lex.ckmeans_exact(x, k, weights=wv)
        d = lex.lex_hybrid(x, k, weights=wv)
        p = lex.lex_plain(x, k, weights=wv)
        tol = 1e-12 * max(1.0, e.sse)
        assert d.sse >= e.sse - tol and p.sse >= e.sse - tol
        assert d.sse <= e.sse * (1 + 1e-4) + tol, (name, k, d.sse, e.sse)
        if len(np.unique(x)) <= lex.DEFAULT_COARSE_GROUPS:
            np.testing.assert_array_equal(d.bounds, e.bounds)  # coarse start = the DP itself


def test_exact_when_atoms_at_most_coarse_groups():
    """atoms <= M: lex_hybrid is ckmeans_exact (same bounds and SSE) and Refine makes 0 moves.
    atoms = M + 1 (one Collapse merge) is only guaranteed SSE >= exact."""
    rng = np.random.default_rng(4096)
    sets = {
        "evenly spaced ties": np.arange(1.0, 9.0),
        "duplicates, d = 5": np.array([0, 0, 0, 1, 1, 2, 2, 2, 2, 5, 9, 9], dtype=float),
        "integers 0..39": rng.integers(0, 40, size=3000).astype(float),
        "gaussian, d = 300": rng.normal(size=300),
        "gaussian, d = 4096 (default M)": np.tile(rng.normal(size=4096), 2),
    }
    for name, x in sets.items():
        d = len(np.unique(x))
        for k in (2, 3, 7, 16, 64):
            if k >= d:
                continue
            e = lex.ckmeans_exact(x, k)
            for m in sorted({d, d + 1, lex.DEFAULT_COARSE_GROUPS}):
                if m < d:
                    continue
                r = lex.lex_hybrid(x, k, coarse_groups=m)
                assert r.start == "coarse_exact", (name, k, m)
                np.testing.assert_array_equal(r.bounds, e.bounds)
                assert r.sse == e.sse and r.moves == 0, (name, k, m, r.sse, e.sse, r.moves)
            c = lex.lex_hybrid(x, k, coarse_groups=d - 1)
            assert c.sse >= e.sse * (1 - 1e-12), (name, k)


def test_ties_follow_r():
    # {0,1,2}, k=2: both splits have SSE 0.5; Ckmeans.1d.dp gives sizes 2|1 (largest cut)
    for f in (lex.lex_hybrid, lex.lex_plain, lex.ckmeans_exact):
        assert list(f([0.0, 1.0, 2.0], 2).sizes) == [2, 1]
    assert list(lex.lex_plain(np.arange(10.0), 3).sizes) == [4, 3, 3]
    assert list(lex.ckmeans_exact(np.arange(10.0), 3).sizes) == [4, 3, 3]


def test_weights_and_api():
    rng = np.random.default_rng(1)
    x = rng.normal(size=3000)
    ones = lex.lex_hybrid(x, 8, weights=np.ones_like(x))
    unweighted = lex.lex_hybrid(x, 8)
    assert ones.sse == pytest.approx(unweighted.sse, rel=1e-12)
    w = rng.lognormal(size=x.size)
    r = lex.ckmeans_exact(x, 8, weights=w)
    lab = r.assign(np.sort(x))
    assert np.array_equal(np.bincount(lab, minlength=8), r.sizes)
    assert r.bounds.dtype == np.int64 and r.centers.dtype == np.float64
    assert lex.lex_hybrid(np.sort(x), 8).presorted and not lex.lex_hybrid(x[::-1], 8).presorted
    # input is not modified; lists and float32 are converted
    y = x.copy(); lex.lex_hybrid(y, 8); assert np.array_equal(x, y)
    assert lex.lex_hybrid(x.astype(np.float32), 8).k == 8
    assert lex.lex_hybrid(list(x[:100]), 4).k == 4


@pytest.mark.parametrize("bad, msg", [
    ([], "empty"), ([1.0, np.nan], "NaN"), ([1.0, np.inf], "Infinity"),
])
def test_errors(bad, msg):
    with pytest.raises(ValueError, match=msg):
        lex.lex_hybrid(np.array(bad, dtype=float), 1)


def test_k_above_atoms_is_clamped():
    x = np.array([3.0, 1.0, 1.0, 2.0, 3.0, 3.0, 2.0])
    w = np.array([1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0])
    for f, kw in ((lex.lex_hybrid, {}), (lex.lex_plain, {}), (METHODS["deck"], {}), (lex.lex_hybrid, {"weights": w})):
        with pytest.warns(UserWarning, match=r"k = 10 exceeds the number of (distinct values|groups at precision 0.0001) \(3\); using k = 3"):
            r = f(x, 10, **kw)
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # no warning when k fits
            r3 = f(x, 3, **kw)
        assert r.k == 3 and list(r.sizes) == [2, 2, 3] and r.sse == 0.0
        assert np.array_equal(r.bounds, r3.bounds) and np.array_equal(r.centers, r3.centers)
    with pytest.warns(UserWarning, match="using k = 3"):
        assert lex.ckmeans_exact(x, 5).k == 3
    with pytest.warns(UserWarning, match="using k = 3"):
        assert lex.ckmeans_exact(x, 5, weights=w).k == 3
    with pytest.warns(UserWarning, match="using k = 1"):
        assert lex.ckmeans_exact([1.0, 1.0], 2).k == 1
    with pytest.warns(UserWarning, match=r"groups at precision 0.1 \(2\)"):
        assert lex.lex_plain([1.0, 1.05, 2.0], 3, precision=0.1).k == 2  # precision groups count as atoms


def test_more_errors():
    with pytest.raises(ValueError, match="weights size"):
        lex.lex_hybrid([1.0, 2.0], 1, weights=[1.0])
    with pytest.raises(ValueError, match="coarse_groups"):
        lex.lex_hybrid([1.0, 2.0], 1, coarse_groups=-1)
    with pytest.raises(TypeError):
        lex.lex_plain([1.0, 2.0], 1, coarse_groups=8)  # lex_plain has no coarse start
    with pytest.raises(ValueError):
        lex.lex_hybrid(np.zeros((2, 2)), 1)


def test_overflow_input_rule_names_the_limit():
    # n * max |x - mean| (with weights max(W, 1) * max |x - mean|, and W) must be <= sqrt(DBL_MAX) / 4
    limit = np.sqrt(np.finfo(np.float64).max) / 4
    assert limit == 3.351951982485649e153
    lim = r"exceeds the limit sqrt\(DBL_MAX\) / 4 = 3\.351951982485649e\+153"
    a = limit / 2  # {-a, a} has mean 0, so n * max |x - mean| = 2a = limit
    up = np.nextafter(a, np.inf)
    for f in (lex.lex_hybrid, lex.lex_plain, lex.ckmeans_exact):
        assert np.isfinite(f([-a, a], 2).sse)
        with pytest.raises(ValueError, match=r"n \* max \|x - mean\| = 3\.3519519824856493e\+153 " + lim
                           + r" \(n = 2, max \|x - mean\| = 1\.6759759912428246e\+153\)"):
            f([-up, up], 2)
    with pytest.raises(ValueError, match=lim + r" \(n = 1000, max \|x - mean\| = 1\.0474849945267654e\+152\)"):
        lex.lex_hybrid(np.repeat([-2.0**505, 2.0**505], 500), 2)
    with pytest.raises(ValueError, match="mean is not finite"):
        lex.lex_hybrid([1e308, 1e308], 1)
    b = limit / 8  # W = 8
    assert np.isfinite(lex.lex_hybrid([-b, b], 2, weights=[4.0, 4.0]).sse)
    with pytest.raises(ValueError, match=r"max\(W, 1\) \* max \|x - mean\| = 3\.3519519824856493e\+153 " + lim):
        lex.lex_hybrid([-np.nextafter(b, np.inf), np.nextafter(b, np.inf)], 2, weights=[4.0, 4.0])
    with pytest.raises(ValueError, match=r"the total weight W = 6\.703903964971298e\+153 " + lim):
        lex.lex_hybrid([0.0, 1.0], 2, weights=[limit, limit])
