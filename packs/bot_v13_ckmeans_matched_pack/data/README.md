# Inputs (large files are not shipped)

## Synthetic data (generated on the Mac)

`bash run/ckmeans_matched.sh` runs `run/gen_synth.py` and checks every needed
file against `data/synth/MANIFEST.sha256`. Do not copy 8M binaries into the
pack; the script writes them under `data/synth/`. Needs NumPy 1.17+ (PCG64 /
`default_rng`). The pinned hashes were produced with NumPy 2.5.5, the same
files as the v10 scaling pack.

Families (same parameters as the lex1d paper/gold `generate.py` and v10):

  uniform    U(0, 1)
  normal     N(0, 1)
  lognormal  LogNormal(0, 1)

This pack uses N = 2^20 and 2^23. Each (family, N) is drawn independently from
`SeedSequence(2048, spawn_key=(stream_id, N))` with stream ids 0, 1, 4. These
are the same draws as the corresponding v10 cells, not prefixes of a longer
file and not the archived Table 4 slices.

## Audio (optional, but expected for the full run)

Place the paper's Dharini Abhogi recording in this `data/` folder, the same
file as v10:

    data/dharini-abhogi.txt.gz
      36,885,053 bytes (gzip -9 -n)
      SHA-256: 3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27

An already-uncompressed `data/dharini-abhogi.txt` is also accepted. The
uncompressed text must have SHA-256

      b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e

N = 8,168,256. If the file is missing, the audio cells are skipped and the
synthetic cells still run.
