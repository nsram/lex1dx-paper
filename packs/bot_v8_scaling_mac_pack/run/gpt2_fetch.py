"""Fetch two unused GPT-2 Large weight tensors (standard library only).

Same source as the paper's GPT-2 section (LEX-9): Hugging Face
openai-community/gpt2-large (MIT), file model.safetensors at the pinned commit
below. Default: HTTPS range reads of the header and the two tensors only
(~46 MB), not the whole 3.2 GB file. --full instead downloads model.safetensors
and checks its SHA-256 against the Hugging Face LFS hash, then extracts.

Tensors (not the paper's h.0.mlp.c_fc.weight):
  h.17.attn.c_attn.weight  float32  1280 x 3840  = 4,915,200 values
  h.35.mlp.c_proj.weight   float32  5120 x 1280  = 6,553,600 values

Writes data/<name>.f32le (little-endian float32, row-major as stored) and
data/gpt2_oos_stats.json. Looks only in this pack's data/ folder.
"""
import hashlib
import json
import math
import os
import struct
import sys
import urllib.request
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(os.path.dirname(HERE), "data")
REPO = "openai-community/gpt2-large"
REVISION = "32b71b12589c2f8d625668d2335a01cac3249519"
URL = f"https://huggingface.co/{REPO}/resolve/{REVISION}/model.safetensors"
FILE_SIZE = 3_247_159_078
FILE_SHA256 = "5f47f3e12f91cd33b662ce7e433b6150ad5512b5884a2cee961b50e9c3bbebce"
UA = "lex1d-paper/gpt2_oos_fetch"
M_MAX = 16384  # default M at k=512; D must be larger so Collapse still runs

TENSORS = [
    dict(name="h.17.attn.c_attn.weight",
         shape=[1280, 3840],
         out="gpt2-large_h17_attn_c_attn_weight.f32le",
         sha256="04c28282e96a10b1375643345351f09dbf3874aaf22a63d87d968fa088b42a72"),
    dict(name="h.35.mlp.c_proj.weight",
         shape=[5120, 1280],
         out="gpt2-large_h35_mlp_c_proj_weight.f32le",
         sha256="283c5aefc44e1c36ec8862e47589dec3fda4465f7c43618aa6dedd614f503836"),
]


def get_range(url, start, end):
    """Inclusive byte range, in chunks."""
    out = bytearray()
    pos = start
    while pos <= end:
        chunk_end = min(end, pos + (1 << 23) - 1)
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Range": f"bytes={pos}-{chunk_end}"})
        with urllib.request.urlopen(req, timeout=180) as r:
            data = r.read()
        if not data:
            sys.exit(f"short read at byte {pos}")
        out.extend(data)
        pos += len(data)
    return bytes(out)


def stats(vals):
    n = len(vals)
    mean = math.fsum(vals) / n
    m2 = math.fsum((x - mean) ** 2 for x in vals) / n
    m4 = math.fsum((x - mean) ** 4 for x in vals) / n
    return dict(min=min(vals), max=max(vals), mean=mean, sd=math.sqrt(m2 * n / (n - 1)),
                kurtosis_excess=m4 / m2 ** 2 - 3.0)


def main():
    full = "--full" in sys.argv[1:]
    os.makedirs(OUTDIR, exist_ok=True)

    # If both tensors are already present and hash-match, just rewrite stats and exit.
    already = True
    for spec in TENSORS:
        path = os.path.join(OUTDIR, spec["out"])
        if not (os.path.isfile(path) and os.path.getsize(path) == 4 * spec["shape"][0] * spec["shape"][1]):
            already = False
            break
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != spec["sha256"]:
            sys.exit(f"{spec['out']} SHA-256 {h.hexdigest()} != pinned {spec['sha256']}")
    # still need stats; fall through to compute them from local files if already

    if already and not full:
        read = None
        file_sha_verified = False
        header_n = None
        header = None
        src_mode = "local-data"
    elif full:
        path = os.path.join(OUTDIR, "model.safetensors")
        if not (os.path.isfile(path) and os.path.getsize(path) == FILE_SIZE):
            print(f"downloading {URL} ({FILE_SIZE/1e9:.2f} GB) to data/model.safetensors", flush=True)
            req = urllib.request.Request(URL, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r, open(path + ".part", "wb") as f:
                while chunk := r.read(1 << 24):
                    f.write(chunk)
            os.replace(path + ".part", path)
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 24), b""):
                h.update(chunk)
        if h.hexdigest() != FILE_SHA256:
            sys.exit(f"model.safetensors SHA-256 {h.hexdigest()} != pinned {FILE_SHA256}")
        print("model.safetensors SHA-256 matches the Hugging Face LFS hash", flush=True)

        def read(a, b):
            with open(path, "rb") as f:
                f.seek(a)
                return f.read(b - a + 1)
        file_sha_verified = True
        src_mode = "full-file"
        header_n, header = None, None
    else:
        read = lambda a, b: get_range(URL, a, b)
        file_sha_verified = False
        src_mode = "range"
        header_n, header = None, None

    recs = []
    if src_mode != "local-data":
        header_n = struct.unpack("<Q", read(0, 7))[0]
        header = json.loads(read(8, 8 + header_n - 1))

    for spec in TENSORS:
        out_path = os.path.join(OUTDIR, spec["out"])
        N = spec["shape"][0] * spec["shape"][1]
        if src_mode == "local-data":
            with open(out_path, "rb") as f:
                raw = f.read()
            tsha = hashlib.sha256(raw).hexdigest()
            offsets = None
        else:
            t = header.get(spec["name"])
            if not t:
                sys.exit(f"tensor {spec['name']} not in the safetensors header")
            if t["dtype"] != "F32" or t["shape"] != spec["shape"]:
                sys.exit(f"{spec['name']}: dtype {t['dtype']} shape {t['shape']}, expected F32 {spec['shape']}")
            a, b = t["data_offsets"]
            if b - a != 4 * N:
                sys.exit(f"{spec['name']}: {b - a} bytes, expected {4 * N}")
            raw = read(8 + header_n + a, 8 + header_n + b - 1)
            tsha = hashlib.sha256(raw).hexdigest()
            if tsha != spec["sha256"]:
                sys.exit(f"{spec['name']} SHA-256 {tsha} != pinned {spec['sha256']}")
            with open(out_path, "wb") as f:
                f.write(raw)
            offsets = [a, b]
        vals32 = array("f")
        vals32.frombytes(raw)
        if sys.byteorder != "little":
            vals32.byteswap()
        vals = array("d", vals32)
        if any(math.isnan(x) or math.isinf(x) for x in vals):
            sys.exit(f"{spec['name']} has NaN or Inf values")
        D = len(set(vals))
        if D <= M_MAX:
            sys.exit(f"{spec['name']}: D = {D} <= {M_MAX}: Collapse at default M would not run")
        recs.append(dict(
            source=f"https://huggingface.co/{REPO}", revision=REVISION,
            file="model.safetensors", file_size=FILE_SIZE, file_sha256=FILE_SHA256,
            file_sha256_verified_here=file_sha_verified, fetch_mode=src_mode,
            tensor=spec["name"], dtype="F32", shape=spec["shape"],
            data_offsets=offsets, order="row-major (as stored)",
            tensor_bytes_sha256=tsha, N=N, D=D, **stats(vals),
            out=spec["out"],
        ))
        print(f"ok {spec['out']} N={N} D={D} sha256={tsha}", flush=True)

    with open(os.path.join(OUTDIR, "gpt2_oos_stats.json"), "w") as f:
        json.dump(recs, f, indent=1)
    print(f"wrote data/gpt2_oos_stats.json and {len(recs)} tensors")


if __name__ == "__main__":
    main()
