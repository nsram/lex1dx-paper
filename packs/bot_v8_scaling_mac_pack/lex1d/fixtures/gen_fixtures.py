"""Generate the shared fixture set (synthetic data only) used by the C++, Python and Scala cross-checks.

    python fixtures/gen_fixtures.py            # checks the committed fixtures/data against SHA256SUMS
    python fixtures/gen_fixtures.py --update   # regenerate fixtures/data, cases.txt and SHA256SUMS

The data files are committed and are the reference. Regenerating is not portable: random.gauss,
lognormvariate and expovariate call the platform's libm (log, cos, exp), whose last-bit rounding
differs between glibc, macOS and Windows. The committed files were generated on Linux (glibc).
Values are written with repr(), which round-trips doubles exactly, with "\n" line endings.
"""
import hashlib
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def bf16(x: float) -> float:
    """Round to bfloat16 precision (8 significant bits), round half to even."""
    if x == 0:
        return 0.0
    m, e = math.frexp(x)
    return math.ldexp(round(m * 256) / 256, e)


def sets():
    r = random.Random(20261115)
    yield "gauss_500", [r.gauss(0, 1) for _ in range(500)]
    yield "uniform_2000", [r.random() for _ in range(2000)]
    yield "lognormal_1000", [r.lognormvariate(0, 1) for _ in range(1000)]
    yield "exponential_800", [r.expovariate(1.0) for _ in range(800)]
    yield "round2_1500", [round(10.2 + 0.15 * r.gauss(0, 1), 2) for _ in range(1500)]
    yield "int30_1000", [float(min(30, int(r.expovariate(0.2)))) for _ in range(1000)]
    yield "grid100_1000", [float(i % 100) for i in range(1000)]
    yield "bf16_1000", [bf16(r.gauss(0, 1)) for _ in range(1000)]
    yield "shifted_1e6_600", [1e6 + r.gauss(0, 1) for _ in range(600)]
    yield "gauss_17000", [r.gauss(0, 1) for _ in range(17000)]  # > 4096 distinct: real coarse step
    yield "micro_3", [0.0, 1.0, 2.0]
    yield "micro_5", [0.0, 1.0, 2.0, 3.0, 4.0]
    yield "micro_6", [2.0, 2.0, 0.0, 1.0, 1.0, 0.0]
    yield "micro_10", [float(i) for i in range(10)]


WEIGHTED = {"gauss_500", "round2_1500", "int30_1000", "gauss_17000"}
ZERO_WEIGHTS = "round2_1500"  # also written with every 5th weight set to 0 (<name>_w0)
# k above the number of distinct values: clamped to it in every language (as Ckmeans.1d.dp)
CLAMPED = [("micro_3", 8, "-"), ("micro_6", 4, "-"), ("grid100_1000", 128, "-"), ("int30_1000", 64, "int30_1000_w")]
KS = [2, 3, 4, 8, 16, 64]
# 60 uniform values (random.random() is exact and portable) whose k = 4 result changes in all four
# presets when the core is built with FMA contraction (x86-64 -ffp-contract=fast -march=haswell, at
# -O2 and at -O3). The other cases did not catch the -O3 build. One case only, to stay small.
FMA = ("fma_60", 29, 60, 4)  # name, seed, n, k


def fma_values():
    r = random.Random(FMA[1])
    return [r.random() for _ in range(FMA[2])]


def write(path, lines):
    with open(path, "w", newline="\n", encoding="ascii") as f:  # "\n" on Windows too
        f.write("\n".join(lines) + "\n")


def checksums():
    sums = []
    for name in sorted(os.listdir(DATA)):
        with open(os.path.join(DATA, name), "rb") as g:
            sums.append(f"{hashlib.sha256(g.read()).hexdigest()}  data/{name}")
    return sums


def read(path):
    with open(path, encoding="ascii") as f:
        return [line.rstrip("\r\n") for line in f]


def check():
    bad = sorted(set(checksums()) ^ set(read(os.path.join(HERE, "SHA256SUMS"))))
    if bad:
        print("fixture check FAILED (fixtures/data differs from SHA256SUMS):", *bad, sep="\n  ")
        return 1
    print(f"ok: {len(read(os.path.join(HERE, 'SHA256SUMS')))} data files match SHA256SUMS; "
          f"{len(read(os.path.join(HERE, 'cases.txt'))) - 1} cases")
    return 0


def update():
    os.makedirs(DATA, exist_ok=True)
    rw = random.Random(7)
    cases = []
    for name, xs in sets():
        write(os.path.join(DATA, name + ".txt"), [repr(float(x)) for x in xs])
        d = len(set(xs))
        for k in KS:
            if k <= d and k < len(xs):
                cases.append(f"{name} {k} -")
        if name in WEIGHTED:
            ws = [rw.lognormvariate(0, 1) for _ in xs]
            write(os.path.join(DATA, name + "_w.txt"), [repr(w) for w in ws])
            for k in (2, 8, 16):
                cases.append(f"{name} {k} {name}_w")
            if name == ZERO_WEIGHTS:
                write(os.path.join(DATA, name + "_w0.txt"), [repr(0.0 if i % 5 == 0 else w) for i, w in enumerate(ws)])
                for k in (2, 8, 16):
                    cases.append(f"{name} {k} {name}_w0")
    cases += [f"{name} {k} {w}" for name, k, w in CLAMPED]
    write(os.path.join(DATA, FMA[0] + ".txt"), [repr(x) for x in fma_values()])
    cases.append(f"{FMA[0]} {FMA[3]} -")
    case_lines = ["# data k weights('-' = none); every case runs presets hybrid, plain, deck, exact"] + cases
    sums = checksums()
    write(os.path.join(HERE, "cases.txt"), case_lines)
    write(os.path.join(HERE, "SHA256SUMS"), sums)
    print(f"updated fixtures/data, cases.txt ({len(cases)} cases) and SHA256SUMS ({len(sums)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(update() if "--update" in sys.argv[1:] else check())
