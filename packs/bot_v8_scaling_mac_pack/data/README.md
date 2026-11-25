# Inputs (large files are not shipped)

## Synthetic data (generated on the Mac)

`bash run/scaling.sh` runs `run/gen_synth.py` and checks every file against
`data/synth/MANIFEST.sha256`. Do not copy 8M binaries into the pack; the
script writes them under `data/synth/`. Needs NumPy 1.17+ (PCG64 /
`default_rng`). The pinned hashes were produced with NumPy 2.5.5.

Families (same parameters as the lex1d paper/gold `generate.py`):

  uniform    U(0, 1)
  normal     N(0, 1)
  lognormal  LogNormal(0, 1)

N = 2^16 .. 2^23, each (family, N) drawn independently from
`SeedSequence(2048, spawn_key=(stream_id, N))` with stream ids 0, 1, 4.

## Audio (optional, but expected for the full run)

Place the paper's Dharini Abhogi recording in this `data/` folder:

    data/dharini-abhogi.txt.gz
      36,885,053 bytes (gzip -9 -n)
      SHA-256: 3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27

An already-uncompressed `data/dharini-abhogi.txt` is also accepted. The
uncompressed text must have SHA-256

      b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e

If the file is missing, the audio cells are skipped and the synthetic (and
GPT-2, if present) cells still run.

## GPT-2 Large tensors (optional; same convention as the v5/v6 packs)

The script will HTTPS range-read the pinned Hugging Face file
(`openai-community/gpt2-large`, commit
`32b71b12589c2f8d625668d2335a01cac3249519`; header + two tensors, ~46 MB,
not the whole 3.2 GB file). Or copy these little-endian float32 files here:

    data/gpt2-large_h17_attn_c_attn_weight.f32le
      SHA-256: 04c28282e96a10b1375643345351f09dbf3874aaf22a63d87d968fa088b42a72

    data/gpt2-large_h35_mlp_c_proj_weight.f32le
      SHA-256: 283c5aefc44e1c36ec8862e47589dec3fda4465f7c43618aa6dedd614f503836

If neither fetch nor copies succeed, h17 and h35 are skipped.
