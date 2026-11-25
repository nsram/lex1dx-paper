#!/usr/bin/env python3
"""Generate small correctness and scale datasets. Deterministic NumPy Generator."""
import os
from pathlib import Path
import numpy as np

OUT = Path("/workspace/lexicon/work/gronlund_bench/data")
OUT.mkdir(parents=True, exist_ok=True)

def write(path, x):
    path = Path(path)
    np.savetxt(path, x, fmt="%.17g")
    print(f"wrote {path} n={len(x)} min={x.min():.6g} max={x.max():.6g} unique={len(np.unique(x))}")

# Correctness: N in {1e3, 1e4}
# Seed scheme: (kind, N) -> stable seed
for n in (1000, 10000):
    write(OUT / f"uniform_n{n}.txt", np.random.Generator(np.random.PCG64(1)).uniform(0.0, 1.0, n))
    write(OUT / f"normal_n{n}.txt", np.random.Generator(np.random.PCG64(2)).normal(0.0, 1.0, n))
    write(OUT / f"lognormal_n{n}.txt", np.random.Generator(np.random.PCG64(3)).lognormal(0.0, 1.0, n))
    # Tie-heavy integers: 20 distinct values, many repeats, not already sorted
    vals = np.random.Generator(np.random.PCG64(4)).integers(0, 20, size=n)
    write(OUT / f"ties_n{n}.txt", vals.astype(float))

# Scale: normal N in {1e5, 1e6}
for n in (100000, 1000000):
    write(OUT / f"normal_n{n}.txt", np.random.Generator(np.random.PCG64(2)).normal(0.0, 1.0, n))

# Audio prefixes from the paper file (read-only original). Not a paper-defined slice.
# First N lines of dharini-abhogi.txt (N=1e5 and 1e6). Full file is 8,168,256 rows.
src = Path("/workspace/lexicon/originals/lex_zip/extracted/lex/empirical/input/dharini-abhogi.txt")
# stream first 1e6 lines only
need = 1000000
chunk = []
with src.open() as f:
    for i, line in enumerate(f):
        if i >= need:
            break
        chunk.append(line)
        if i + 1 in (100000, 1000000):
            dest = OUT / f"audio_prefix_n{i+1}.txt"
            dest.write_text("".join(chunk))
            # count unique later if cheap
            print(f"wrote {dest} n={i+1} (prefix of {src})")

# unique counts for audio prefixes
for n in (100000, 1000000):
    p = OUT / f"audio_prefix_n{n}.txt"
    x = np.loadtxt(p)
    print(f"audio_prefix_n{n}: unique={len(np.unique(x))} min={x.min():.6g} max={x.max():.6g}")
