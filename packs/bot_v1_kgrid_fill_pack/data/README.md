# Inputs (large files are not shipped)

## Synthetic data (generated on the Mac)

`bash run/kgrid_fill.sh` runs `run/gen_synth.py` and checks every full-run
file against `data/synth/MANIFEST.sha256`. Do not copy 8M binaries into
the pack; the script writes them under `data/synth/`. Needs NumPy 1.17+
(PCG64 / `default_rng`). The pinned hashes were produced with NumPy 2.5.5.

The full run generates only N = 2^20 and 2^23. A `--smoke` run also
generates N = 4096, which is not in the shipped manifest.

Families (same parameters as the paper's generator):

  uniform    U(0, 1)
  normal     N(0, 1)
  lognormal  LogNormal(0, 1)

Each (family, N) is drawn independently from
`SeedSequence(2048, spawn_key=(stream_id, N))` with stream ids 0, 1, 4.

## Audio (optional, but expected for the full run)

Place the Dharini Abhogi recording in this `data/` folder:

    data/dharini-abhogi.txt.gz
      36,885,053 bytes (gzip -9 -n)
      SHA-256: 3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27

An already-uncompressed `data/dharini-abhogi.txt` is also accepted. The
uncompressed text must have SHA-256

      b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e

If the file is missing, the audio cells are skipped and the synthetic
cells still run.
