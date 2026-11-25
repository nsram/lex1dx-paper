#!/usr/bin/env python3
"""Generate synthetic 1-D float64 datasets matching the lex1d paper/gold families.

Families and parameters (from lex1d-paper/data/generate.py):

  uniform    U(0, 1)                  paper: uniform_continuous, stream 0
  normal     N(0, 1)                  paper: gaussian_continuous, stream 1
  lognormal  LogNormal(mu=0, sigma=1) paper: lognormal_continuous, stream 4

Default N grid: 2^16 .. 2^23 (every power, 8 levels). Each (family, N) is drawn
independently (a smaller N is not a prefix of a larger file).

Seeding (deterministic across machines on the same NumPy Generator/PCG64):

  numpy.random.SeedSequence(2048, spawn_key=(stream_id, N))
  rng = numpy.random.default_rng(ss)   # PCG64; 2048 is the paper generate.py seed

Requires NumPy 1.17+ (Generator, PCG64, default_rng). Bit-identical files are
expected for a given NumPy major.minor; this box used 2.5.5. Do not use the
legacy RandomState/MT19937.

Writes little-endian IEEE-754 float64 binaries:

  data/synth/<family>_n<N>.f64

and a manifest data/synth/MANIFEST.sha256 (plus MANIFEST.json).

Usage (from the radix_sort folder, or any cwd with --out):

  python3 gen_synth.py
  python3 gen_synth.py --exp 16-20
  python3 gen_synth.py --exp 16,18,23 --out data/synth
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

try:
    import numpy as np
except ImportError as e:
    sys.exit("gen_synth.py needs NumPy 1.17+ (PCG64 / default_rng): " + str(e))

# Paper generate.py
PAPER_SEED = 2048
STREAM_ID = {"uniform": 0, "normal": 1, "lognormal": 4}
FAMILIES = ("uniform", "normal", "lognormal")
DEFAULT_EXP_LO = 16
DEFAULT_EXP_HI = 23  # inclusive; 8 levels


def parse_exp(s: str) -> list[int]:
    s = s.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        lo, hi = int(a), int(b)
        if lo > hi:
            lo, hi = hi, lo
        return list(range(lo, hi + 1))
    return [int(x) for x in s.split(",") if x.strip() != ""]


def make_rng(family: str, n: int) -> np.random.Generator:
    if family not in STREAM_ID:
        raise ValueError(f"unknown family {family}")
    ss = np.random.SeedSequence(PAPER_SEED, spawn_key=(STREAM_ID[family], int(n)))
    return np.random.default_rng(ss)  # PCG64


def draw(family: str, n: int) -> np.ndarray:
    rng = make_rng(family, n)
    if family == "uniform":
        x = rng.uniform(0.0, 1.0, n)
    elif family == "normal":
        x = rng.normal(0.0, 1.0, n)
    elif family == "lognormal":
        x = rng.lognormal(0.0, 1.0, n)
    else:
        raise ValueError(family)
    if x.dtype != np.float64:
        x = x.astype(np.float64, copy=False)
    return np.ascontiguousarray(x, dtype=np.float64)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if tuple(int(p) for p in np.__version__.split(".")[:2]) < (1, 17):
        sys.exit(f"NumPy {np.__version__} is too old; need 1.17+ for default_rng/PCG64")

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exp", default=f"{DEFAULT_EXP_LO}-{DEFAULT_EXP_HI}",
                   help="exponents for N=2^e, inclusive range or comma list (default 16-23)")
    p.add_argument("--out", default=None, help="output directory (default: <script>/data/synth or ./data/synth)")
    p.add_argument("--families", default="uniform,normal,lognormal",
                   help="comma list of families (default: all three)")
    p.add_argument("--no-write-manifest", action="store_true",
                   help="generate files only; do not overwrite a shipped MANIFEST")
    args = p.parse_args()

    exponents = parse_exp(args.exp)
    if not exponents or min(exponents) < 0 or max(exponents) > 32:
        sys.exit("bad --exp")
    families = tuple(x.strip() for x in args.families.split(",") if x.strip())
    for f in families:
        if f not in STREAM_ID:
            sys.exit(f"unknown family {f}; expected one of {FAMILIES}")

    if args.out:
        out = Path(args.out)
    else:
        here = Path(__file__).resolve().parent
        cand = here / "data" / "synth"
        out = cand if (here / "data").is_dir() or here.name == "radix_sort" else Path("data/synth")
    out.mkdir(parents=True, exist_ok=True)

    recs = []
    for family in families:
        for e in exponents:
            n = 1 << e
            path = out / f"{family}_n{n}.f64"
            x = draw(family, n)
            x = np.ascontiguousarray(x, dtype="<f8")
            path.write_bytes(x.tobytes())
            digest = sha256_file(path)
            rec = dict(
                family=family,
                N=n,
                exponent=e,
                path=str(path.name),
                bytes=n * 8,
                sha256=digest,
                seed_root=PAPER_SEED,
                spawn_key=[STREAM_ID[family], n],
                rng="numpy.random.default_rng(SeedSequence(2048, spawn_key=(stream_id, N)))  # PCG64",
                params={"uniform": "U(0,1)", "normal": "N(0,1)", "lognormal": "LogNormal(0,1)"}[family],
                numpy=np.__version__,
            )
            recs.append(rec)
            print(f"ok {path.name} N={n} sha256={digest}", flush=True)

    if args.no_write_manifest:
        print(f"generated {len(recs)} files under {out} (manifest not rewritten)")
        return

    man = out / "MANIFEST.sha256"
    with man.open("w") as f:
        f.write(f"# gen_synth.py  numpy={np.__version__}  seed={PAPER_SEED}  "
                f"spawn_key=(paper_stream_id, N)  rng=default_rng/PCG64\n")
        f.write("# family  N  sha256  file\n")
        for r in recs:
            f.write(f"{r['family']:10s} {r['N']:8d}  {r['sha256']}  {r['path']}\n")
    (out / "MANIFEST.json").write_text(json.dumps(dict(
        numpy=np.__version__,
        seed_root=PAPER_SEED,
        stream_id=STREAM_ID,
        rng="default_rng/PCG64",
        requirement="NumPy 1.17+",
        files=recs,
    ), indent=2) + "\n")
    print(f"wrote {man} and MANIFEST.json ({len(recs)} files)")


if __name__ == "__main__":
    main()
