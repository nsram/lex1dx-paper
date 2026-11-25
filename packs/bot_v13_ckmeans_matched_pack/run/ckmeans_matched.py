#!/usr/bin/env python3
"""Orchestrate the v11 Ckmeans.1d.dp matched rerun. Invoked by run/ckmeans_matched.sh."""
from __future__ import annotations

import csv
import hashlib
import os
import re
import select
import shutil
import subprocess
import sys
import time
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
DATA = PKG / "data"
WORK = PKG / "work"
OUT = PKG / "out" / "ckmeans_matched"
BUILD = PKG / "build"
REF = PKG / "ref" / "v10_wilber_sse.csv"

# Modest grid, matched to v10 inputs. N = 2^20 and 2^23 as in the v10 scaling
# grid; k covers v10's {64,256,512} plus 16 and 128.
KS = (16, 64, 128, 256, 512)
FAMILIES = ("uniform", "normal", "lognormal")
FULL_N = (1_048_576, 8_388_608)
SMOKE_N = (4_096,)
SMOKE_KS = (4, 16)

AUDIO_GZ_NAME = "dharini-abhogi.txt.gz"
AUDIO_GZ_SHA = "3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27"
AUDIO_TXT_SHA = "b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e"
AUDIO_N = 8_168_256

# Peak-RSS estimator for Ckmeans.1d.dp method=linear, fit to the archived
# 4.3.6 runs at N=2^20 (k=64 ~1.20 GiB, k=256 ~4.22 GiB, k=512 ~8.25 GiB):
#   RSS ≈ 0.22 GiB + 16 bytes * N * k
# Time estimator: archived sorted-fit median ~21.6 s at N=2^20, k=512, O(kn).
BYTES_PER_NK = 16.0
BASE_GIB = 0.22
SEC_PER_NK = 21.6 / (1_048_576 * 512)

DEFAULT_MEM_CAP_GIB = 24.0
DEFAULT_TIME_CAP_S = 1800.0

CSV_FIELDS = [
    "dataset", "method", "ckmeans_method", "version", "k", "n", "D",
    "input", "status",
    "sort_s", "post_s", "time_min_s", "time_max_s",
    "sse", "tot_withinss", "exact_sse", "ree_pct", "exact_sse_source",
    "k_eff", "warm", "reps",
    "peak_rss_mb", "est_rss_gib", "est_fit_s", "elapsed_s",
    "note",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_cell(text: str) -> dict | None:
    line = next((ln for ln in text.splitlines() if ln.startswith("CELL ")), None)
    if not line:
        return None
    kv = {}
    for p in line.split()[1:]:
        if "=" in p:
            a, b = p.split("=", 1)
            kv[a] = b
    return kv


def load_csv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def append_row(path: Path, row: dict) -> None:
    new = not path.is_file()
    with path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in CSV_FIELDS})
        f.flush()


def done(rows: list[dict], dataset: str, method: str, k: int, input_kind: str) -> bool:
    for r in rows:
        if (r.get("dataset") == dataset and r.get("method") == method
                and int(r.get("k") or 0) == k and (r.get("input") or "") == input_kind):
            if r.get("status") in ("ok", "skip_k_ge_D"):
                return True
    return False


def mem_cap_gib() -> float:
    raw = os.environ.get("LEX_MEM_CAP_GIB", str(DEFAULT_MEM_CAP_GIB))
    try:
        v = float(raw)
    except ValueError:
        sys.exit(f"error: LEX_MEM_CAP_GIB={raw!r} is not a number")
    if v <= 0:
        sys.exit("error: LEX_MEM_CAP_GIB must be positive")
    return v


def time_cap_s() -> float:
    raw = os.environ.get("LEX_TIME_CAP_S", str(DEFAULT_TIME_CAP_S))
    try:
        v = float(raw)
    except ValueError:
        sys.exit(f"error: LEX_TIME_CAP_S={raw!r} is not a number")
    if v <= 0:
        sys.exit("error: LEX_TIME_CAP_S must be positive")
    return v


def estimate_rss_gib(n: int, k: int) -> float:
    return BASE_GIB + (BYTES_PER_NK * n * k) / (1024.0 ** 3)


def estimate_fit_s(n: int, k: int) -> float:
    return SEC_PER_NK * n * k


def find_bin(name: str) -> Path:
    for p in (BUILD / name, BUILD / "Release" / name):
        if p.is_file() and os.access(p, os.X_OK):
            return p
    sys.exit(f"error: missing binary {name} under {BUILD}")


def ru_maxrss_mb(ru) -> float:
    # wait4 / getrusage: bytes on macOS, KiB on Linux (same as v10 scale_cell).
    if sys.platform == "darwin":
        return float(ru.ru_maxrss) / (1024.0 * 1024.0)
    return float(ru.ru_maxrss) / 1024.0


def run_wait4(cmd: list[str], timeout: float | None, env: dict | None = None
              ) -> tuple[int, str, str, float | None, bool]:
    """Run cmd, return (rc, stdout, stderr, peak_rss_mb, timed_out).

    Peak RSS is the child's ru_maxrss from wait4, so it is the cell process
    (the same quantity v10 reports via getrusage).
    """
    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env,
    )
    assert p.stdout is not None and p.stderr is not None
    out = b""
    err = b""
    deadline = (time.time() + timeout) if timeout else None
    timed_out = False
    rusage = None
    status = None
    try:
        while True:
            remaining = None if deadline is None else max(0.0, deadline - time.time())
            if deadline is not None and remaining == 0.0 and p.poll() is None:
                timed_out = True
                p.kill()
            wait_flags = os.WNOHANG
            pid, status, ru = os.wait4(p.pid, wait_flags)
            if pid != 0:
                rusage = ru
                # Drain anything still in the pipes.
                for stream in (p.stdout, p.stderr):
                    try:
                        chunk = os.read(stream.fileno(), 1 << 20)
                    except OSError:
                        chunk = b""
                    if stream is p.stdout:
                        out += chunk
                    else:
                        err += chunk
                break
            to = 0.1 if remaining is None else min(0.1, remaining if remaining > 0 else 0.01)
            ready, _, _ = select.select([p.stdout, p.stderr], [], [], to)
            if p.stdout in ready:
                chunk = os.read(p.stdout.fileno(), 1 << 16)
                if chunk:
                    out += chunk
            if p.stderr in ready:
                chunk = os.read(p.stderr.fileno(), 1 << 16)
                if chunk:
                    err += chunk
    finally:
        if status is None:
            try:
                pid, status, ru = os.wait4(p.pid, 0)
                rusage = ru
            except ChildProcessError:
                status = 0
                if p.returncode is not None:
                    status = p.returncode
    # Decode
    so = out.decode("utf-8", errors="replace")
    se = err.decode("utf-8", errors="replace")
    if os.WIFEXITED(status):
        rc = os.WEXITSTATUS(status)
    elif os.WIFSIGNALED(status):
        rc = 128 + os.WTERMSIG(status)
    else:
        rc = p.returncode if p.returncode is not None else 1
    if timed_out:
        rc = 124
        se = se + ("" if se.endswith("\n") or se == "" else "\n") + f"timeout after {timeout:.0f}s"
    peak = ru_maxrss_mb(rusage) if rusage is not None else None
    return rc, so, se, peak, timed_out


def parse_manifest(path: Path) -> dict[str, str]:
    rec = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            rec[parts[3]] = parts[2]
        elif len(parts) == 2:
            rec[parts[1] if "/" in parts[1] or parts[1].endswith(".f64") else parts[0]] = parts[0]
    return rec


def parse_exp_from_ns(ns: tuple[int, ...]) -> str:
    # gen_synth.py --exp wants exponents, or a range.
    exps = []
    for n in ns:
        e = n.bit_length() - 1
        if (1 << e) != n:
            raise ValueError(f"N={n} is not a power of two")
        exps.append(e)
    return ",".join(str(e) for e in exps)


def generate_synth(python: str, ns: tuple[int, ...], check_manifest: bool) -> None:
    out = DATA / "synth"
    out.mkdir(parents=True, exist_ok=True)
    exp = parse_exp_from_ns(ns)
    expect = parse_manifest(out / "MANIFEST.sha256") if (out / "MANIFEST.sha256").is_file() else {}
    need = []
    for n in ns:
        for fam in FAMILIES:
            fname = f"{fam}_n{n}.f64"
            path = out / fname
            digest = expect.get(fname)
            if path.is_file() and digest and sha256_file(path) == digest:
                continue
            if path.is_file() and not digest:
                continue  # smoke sizes are not in the shipped manifest
            need.append(fname)
    man_path = out / "MANIFEST.sha256"
    json_path = out / "MANIFEST.json"
    man_bytes = man_path.read_bytes() if man_path.is_file() else None
    json_bytes = json_path.read_bytes() if json_path.is_file() else None
    if need or any(not (out / f"{fam}_n{n}.f64").is_file() for fam in FAMILIES for n in ns):
        log(f"generating synthetic data (--exp {exp}) with {python} ...")
        r = subprocess.run(
            [python, str(PKG / "run" / "gen_synth.py"),
             "--exp", exp, "--out", str(out), "--no-write-manifest"],
            check=False,
        )
        if man_bytes is not None:
            man_path.write_bytes(man_bytes)
        if json_bytes is not None:
            json_path.write_bytes(json_bytes)
        if r.returncode != 0:
            sys.exit("error: gen_synth.py failed (need NumPy 1.17+, PCG64/default_rng)")
    else:
        if man_bytes is not None:
            man_path.write_bytes(man_bytes)
        if json_bytes is not None:
            json_path.write_bytes(json_bytes)
    if not check_manifest:
        log(f"synth ready (smoke sizes {list(ns)}, hashes not checked against the v10 manifest)")
        return
    allowed = set(ns)
    bad = []
    for fname, digest in expect.items():
        m = re.search(r"_n(\d+)\.f64$", fname)
        if not m or int(m.group(1)) not in allowed:
            continue
        path = out / fname
        if not path.is_file():
            bad.append(f"{fname} missing")
            continue
        got = sha256_file(path)
        if got != digest:
            bad.append(f"{fname} sha256 {got} != pinned {digest}")
    if bad:
        local = out / "MANIFEST_local.sha256"
        with local.open("w") as f:
            for fname in sorted(expect):
                m = re.search(r"_n(\d+)\.f64$", fname)
                if m and int(m.group(1)) in allowed and (out / fname).is_file():
                    f.write(f"{sha256_file(out / fname)}  {fname}\n")
        log("warning: %d synthetic files differ from the box hashes "
            "(platform libm differences in exp/log); local hashes in %s:\n  %s"
            % (len(bad), local.name, "\n  ".join(bad)))
        return
    log(f"synth ok ({len(allowed)} sizes x 3 families, hashes match MANIFEST.sha256)")


def prepare_audio(txt64: Path) -> Path | None:
    gz = DATA / AUDIO_GZ_NAME
    txt = DATA / "dharini-abhogi.txt"
    work_txt = WORK / "dharini-abhogi.txt"
    audio_f64 = WORK / "audio.f64"
    if audio_f64.is_file() and audio_f64.stat().st_size == AUDIO_N * 8:
        return audio_f64
    if not (gz.is_file() or txt.is_file()):
        log("audio not in data/; skipping (place dharini-abhogi.txt.gz as usual to include it)")
        return None
    if gz.is_file():
        subprocess.run(["gunzip", "-c", str(gz)], check=True, stdout=work_txt.open("wb"))
    else:
        shutil.copy(txt, work_txt)
    got = sha256_file(work_txt)
    if got != AUDIO_TXT_SHA:
        log(f"warning: audio text SHA-256 {got} != {AUDIO_TXT_SHA}; skipping audio")
        return None
    subprocess.run([str(txt64), str(work_txt), str(audio_f64)], check=True)
    if audio_f64.stat().st_size != AUDIO_N * 8:
        sys.exit("error: audio f64 has the wrong size")
    log("audio ready")
    return audio_f64


def load_v10_sse() -> dict[tuple[str, int], str]:
    rec = {}
    if not REF.is_file():
        return rec
    for r in load_csv(REF):
        if r.get("status") == "ok" and r.get("sse"):
            rec[(r["dataset"], int(r["k"]))] = r["sse"]
    return rec


def score_bounds(score: Path, sorted_path: Path, bounds: Path) -> str:
    rc, so, se, _, _ = run_wait4([str(score), "bounds", str(sorted_path), str(bounds)], timeout=120)
    if rc != 0:
        raise RuntimeError(f"score_sse bounds failed rc={rc}: {se or so}")
    m = re.search(r"SSE=([^\s]+)", so)
    if not m:
        raise RuntimeError(f"score_sse bounds printed no SSE: {so}")
    return m.group(1)


def wilber_sse(score: Path, sorted_path: Path, k: int, cache: dict[tuple[str, int], tuple[str, str]],
               dataset: str, log_path: Path) -> tuple[str, str]:
    key = (dataset, k)
    if key in cache and cache[key][0]:
        return cache[key]
    # Prefer a just-computed value. Fall back to shipped v10 SSE* after we try.
    to = max(60.0, 120.0 * (sorted_path.stat().st_size / 8) / 8_000_000)
    log(f"======= wilber SSE* {dataset} k={k} =======")
    t0 = time.time()
    rc, so, se, peak, timed_out = run_wait4([str(score), "wilber", str(sorted_path), str(k)], timeout=to)
    elapsed = time.time() - t0
    with log_path.open("a") as f:
        f.write(f"CMD score_sse wilber {sorted_path} {k} rc={rc} elapsed={elapsed:.3f}\n")
        f.write(so)
        if se:
            f.write(se if se.endswith("\n") else se + "\n")
    if rc == 0:
        m = re.search(r"SSE=([^\s]+)", so)
        if m:
            src = "this pack, score_sse wilber (v10 patched wilber-interp + raw_sse)"
            cache[key] = (m.group(1), src)
            log(f"wilber SSE* {dataset} k={k} {m.group(1)} ({elapsed:.1f}s, rss={peak})")
            return cache[key]
    # Fall back to the shipped v10 value if present (k=64,256,512 only).
    v10 = load_v10_sse()
    if key in v10:
        log(f"warning: wilber SSE* failed for {dataset} k={k}; using shipped v10 value")
        cache[key] = (v10[key], "v10 Mac run (ref/v10_wilber_sse.csv); this-pack wilber failed")
        return cache[key]
    log(f"warning: no wilber SSE* for {dataset} k={k} (rc={rc})")
    cache[key] = ("", "")
    return cache[key]


def write_summary(rows: list[dict], path: Path, smoke: bool,
                  mem_cap: float, time_cap: float, ckmeans_method: str) -> None:
    def fnum(r, k, nd=3):
        try:
            return f"{float(r.get(k) or 0):.{nd}f}"
        except ValueError:
            return r.get(k) or ""

    lines = []
    lines.append("# Ckmeans.1d.dp matched to the lex1dx v10 protocol")
    lines.append("")
    if smoke:
        lines.append(
            "SMOKE: synthetic N = 4096, three families, k ∈ {4,16}, one warm-up and "
            "one timed rep, both pre-sorted and unsorted inputs. Times are that single "
            "rep, not a median. SSE is the v10 two-pass long-double scorer; REE is vs "
            "wilber-interp SSE* computed in this pack with the same scorer."
        )
    else:
        lines.append(
            "One warm-up and 3 timed reps; times are medians, in seconds. "
            "Inputs, seeds, and files are those of the v10 scaling pack "
            "(`run/gen_synth.py`, SeedSequence(2048, spawn_key=(stream_id, N)); "
            "audio is the v10 Dharini Abhogi recording, N = 8,168,256). "
            f"`Ckmeans.1d.dp` ran with `method = \"{ckmeans_method}\"` (fastest exact method, "
            "O(kn) after sort, O(kn) memory). "
            "**Post-sort time (`input=sorted`)** is a call on an array already sorted with "
            "`std::sort` (the same sort as v10, timed separately as `sort_s`). "
            "`Ckmeans.1d.dp` still sorts internally; that internal sort is on already-sorted "
            "data and is **not excluded**. Unsorted fits (`input=raw`) include the internal "
            "sort of the unsorted array, matching the archived Table 4 column. "
            "SSE is scored on the raw sorted values with the same two-pass long-double "
            "`raw_sse` as v10 `scale_cell.cpp`. Exact SSE* is unweighted Grønlund "
            "wilber-interp (interpolation search, the v10 patched solver) scored with the "
            "same function. REE is (SSE − SSE*) / SSE* × 100. "
            "Peak RSS is the cell process (`wait4` / `getrusage`), in MiB. "
            f"A cell is skipped, and recorded as `status=skipped`, if the estimated peak "
            f"RSS exceeds {mem_cap:g} GiB (`LEX_MEM_CAP_GIB`) or the estimated wall time "
            f"for the warm-up plus reps exceeds {time_cap:g} s (`LEX_TIME_CAP_S`). "
            "The estimator is a fit to the archived Ckmeans.1d.dp 4.3.6 runs at N = 2^20: "
            "RSS ≈ 0.22 GiB + 16 bytes × N × k, and time ≈ 21.6 s × (N / 2^20) × (k / 512) "
            "per fit."
        )
    lines.append("")

    sort_s = {}
    by = {}
    for r in rows:
        if r.get("method") == "sort":
            sort_s[r["dataset"]] = r
            continue
        by[(r.get("dataset"), r.get("input"), int(r.get("k") or 0))] = r

    if sort_s:
        lines.append("## Sort reference (`std::sort`, not in post-sort)")
        lines.append("")
        lines.append("| dataset | n | D | sort s | peak RSS MiB |")
        lines.append("|---|---:|---:|---:|---:|")
        for ds, r in sorted(sort_s.items(), key=lambda kv: kv[0]):
            lines.append(f"| {ds} | {r.get('n','')} | {r.get('D','')} | {fnum(r,'sort_s',3)} | {fnum(r,'peak_rss_mb',1)} |")
        lines.append("")

    lines.append("## Ckmeans.1d.dp")
    lines.append("")
    lines.append(
        "| dataset | k | input | status | n | D | post-sort / fit s | "
        "SSE | SSE* | REE % | tot.withinss | RSS MiB | est. RSS GiB | version | note |"
    )
    lines.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")

    def ds_key(d: str):
        if "_n" in d:
            try:
                return (0, int(d.split("_n")[-1]), d)
            except ValueError:
                return (1, 0, d)
        return (1, 0, d)

    datasets = sorted({k[0] for k in by} | set(sort_s), key=ds_key)
    ks = sorted({k[2] for k in by})
    for ds in datasets:
        for k in ks:
            for inp in ("sorted", "raw"):
                r = by.get((ds, inp, k))
                if not r:
                    continue
                ree = r.get("ree_pct") or ""
                try:
                    ree = f"{float(ree):.5g}" if ree not in ("", "nan", "NAN") else ""
                except ValueError:
                    pass
                note = (r.get("note") or "").replace("|", "/")
                lines.append(
                    f"| {ds} | {k} | {inp} | {r.get('status')} | {r.get('n')} | {r.get('D')} | "
                    f"{fnum(r,'post_s',3)} | {r.get('sse') or ''} | {r.get('exact_sse') or ''} | "
                    f"{ree} | {r.get('tot_withinss') or ''} | {fnum(r,'peak_rss_mb',1)} | "
                    f"{fnum(r,'est_rss_gib',2)} | {r.get('version') or ''} | {note} |"
                )
    lines.append("")
    lines.append(
        "See `out/ckmeans_matched/runs.csv` and `out/ckmeans_matched/run.log`. "
        "Re-run the same command to resume; finished (`ok`, `skip_k_ge_D`) cells are skipped; `skipped` cells are re-checked against the current caps."
    )
    path.write_text("\n".join(lines) + "\n")


def r_env(base: dict[str, str], rlib: Path) -> dict[str, str]:
    env = dict(base)
    existing = env.get("R_LIBS_USER", "")
    env["R_LIBS_USER"] = str(rlib) if not existing else str(rlib) + os.pathsep + existing
    # Also prepend to R_LIBS so library() finds the pack-local install first.
    existing2 = env.get("R_LIBS", "")
    env["R_LIBS"] = str(rlib) if not existing2 else str(rlib) + os.pathsep + existing2
    return env


def main() -> None:
    args = sys.argv[1:]
    smoke = "--smoke" in args
    sorted_only = "--sorted-only" in args
    unknown = [a for a in args if a not in ("--smoke", "--sorted-only")]
    if unknown:
        sys.exit("usage: ckmeans_matched.py [--smoke] [--sorted-only]")

    python = os.environ.get("PYTHON_BIN") or sys.executable
    rscript = os.environ.get("RSCRIPT_BIN") or "Rscript"
    ckmeans_method = os.environ.get("LEX_CKMEANS_METHOD") or "linear"
    if ckmeans_method not in ("linear", "loglinear", "quadratic"):
        sys.exit("error: LEX_CKMEANS_METHOD must be linear, loglinear, or quadratic")

    sort_bin = find_bin("sort_f64")
    score = find_bin("score_sse")
    txt64 = find_bin("txt_to_f64")
    cell_r = PKG / "run" / "ckmeans_cell.R"
    if not cell_r.is_file():
        sys.exit(f"error: missing {cell_r}")
    rlib = PKG / "rlib"
    rlib.mkdir(parents=True, exist_ok=True)
    env = r_env(os.environ.copy(), rlib)

    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "runs.csv"
    log_path = OUT / "run.log"
    rows = load_csv(csv_path)

    ns = SMOKE_N if smoke else FULL_N
    ks = SMOKE_KS if smoke else KS
    warm = 1
    reps = 1 if smoke else 3
    mem_cap = mem_cap_gib()
    time_cap = time_cap_s()
    inputs = ("sorted",) if sorted_only else ("sorted", "raw")

    generate_synth(python, ns, check_manifest=not smoke)

    datasets: list[tuple[str, Path]] = []
    for n in ns:
        for fam in FAMILIES:
            p = DATA / "synth" / f"{fam}_n{n}.f64"
            if p.is_file():
                datasets.append((f"{fam}_n{n}", p))
    if not smoke:
        audio = prepare_audio(txt64)
        if audio is not None:
            datasets.append(("audio", audio))
    if not datasets:
        sys.exit("error: no datasets to run")

    log(f"datasets={[d for d,_ in datasets]} ks={list(ks)} reps={reps} warm={warm} "
        f"inputs={list(inputs)} method={ckmeans_method} "
        f"mem_cap_gib={mem_cap:g} time_cap_s={time_cap:g} smoke={int(smoke)}")
    t_all = time.time()
    v10_sse = load_v10_sse()
    sse_cache: dict[tuple[str, int], tuple[str, str]] = {
        k: (v, "v10 Mac run (ref/v10_wilber_sse.csv)") for k, v in v10_sse.items()
    }
    # Always recompute SSE* in this pack when missing from v10 (k=16,128, and smoke).
    # For k present in v10, still recompute so REE uses this machine's input + the
    # same scorer; the v10 value is a cross-check, logged if it differs.

    for name, raw in datasets:
        n = raw.stat().st_size // 8
        sorted_path = WORK / f"{name}.sorted.f64"
        if not done(rows, name, "sort", 0, "") or not sorted_path.is_file():
            log(f"======= sort {name} n={n} =======")
            to = max(30.0, 60.0 * (n / 8_000_000))
            rc, out, err, peak, _ = run_wait4(
                [str(sort_bin), str(raw), str(sorted_path), name, str(reps)],
                timeout=to,
            )
            with log_path.open("a") as f:
                f.write(out + err + f"\nrc={rc}\n")
            kv = parse_cell(out)
            if rc != 0 or not kv:
                row = dict(dataset=name, method="sort", k=0, n=n, input="",
                           status="timeout" if rc == 124 else "crash",
                           note=(err or out)[-400:], reps=reps,
                           peak_rss_mb="" if peak is None else f"{peak:.1f}")
                append_row(csv_path, row)
                rows.append(row)
                log(f"SORT FAILED {name} rc={rc}")
                continue
            row = {k: kv.get(k, "") for k in CSV_FIELDS}
            row["dataset"] = name
            row["method"] = "sort"
            row["n"] = kv.get("n", str(n))
            if peak is not None and not row.get("peak_rss_mb"):
                row["peak_rss_mb"] = f"{peak:.1f}"
            append_row(csv_path, row)
            rows.append(row)
            log(out.strip() if out.strip().startswith("CELL") else next(
                (ln for ln in out.splitlines() if ln.startswith("CELL ")), out))
        else:
            log(f"skip sort {name} (already done)")

        if not sorted_path.is_file():
            continue

        sort_row = next((r for r in rows if r.get("dataset") == name and r.get("method") == "sort"), None)
        sort_s = (sort_row or {}).get("sort_s", "")
        D = (sort_row or {}).get("D", "")

        for k in ks:
            # Compute SSE* once per (dataset, k), on the sorted file.
            need_sse = any(
                not done(rows, name, "Ckmeans.1d.dp", k, inp) for inp in inputs
            )
            exact, exact_src = ("", "")
            if need_sse:
                # Recompute unless we already ran wilber in this process for this key
                # *and* the source is this pack. If only a v10 shipped value exists
                # (k=64/256/512), recompute so this machine's input is the reference.
                have = sse_cache.get((name, k))
                if have and have[0] and have[1].startswith("this pack"):
                    exact, exact_src = have
                else:
                    exact, exact_src = wilber_sse(score, sorted_path, k, sse_cache, name, log_path)
                    # Cross-check against shipped v10 SSE* when both exist.
                    v10 = v10_sse.get((name, k))
                    if v10 and exact and v10 != exact:
                        log(f"note: wilber SSE* for {name} k={k} is {exact}; "
                            f"v10 Mac value was {v10} (inputs or long-double may differ)")

            for inp in inputs:
                if done(rows, name, "Ckmeans.1d.dp", k, inp):
                    log(f"skip {name} Ckmeans k={k} input={inp}")
                    continue
                est_rss = estimate_rss_gib(n, k)
                est_fit = estimate_fit_s(n, k)
                est_wall = (warm + reps) * est_fit + 8.0
                skip_reason = ""
                if est_rss > mem_cap:
                    skip_reason = (
                        f"skipped: estimated peak RSS {est_rss:.1f} GiB > cap {mem_cap:g} GiB "
                        f"(method={ckmeans_method})"
                    )
                elif est_wall > time_cap:
                    skip_reason = (
                        f"skipped: estimated wall time {est_wall:.0f} s > cap {time_cap:g} s "
                        f"(method={ckmeans_method}, {warm}+{reps} fits at ~{est_fit:.1f} s each)"
                    )
                if skip_reason:
                    row = dict(
                        dataset=name, method="Ckmeans.1d.dp", ckmeans_method=ckmeans_method,
                        k=k, n=n, D=D, input=inp, status="skipped",
                        sort_s=sort_s, exact_sse=exact, exact_sse_source=exact_src,
                        warm=warm, reps=reps,
                        est_rss_gib=f"{est_rss:.3f}", est_fit_s=f"{est_fit:.3f}",
                        note=skip_reason,
                    )
                    append_row(csv_path, row)
                    rows.append(row)
                    log(f"SKIP {name} k={k} input={inp}  {skip_reason}")
                    continue

                in_path = sorted_path if inp == "sorted" else raw
                bounds = WORK / f"{name}_k{k}_{inp}.bounds"
                if bounds.exists():
                    bounds.unlink()
                cmd = [
                    rscript, str(cell_r), str(in_path), str(k), ckmeans_method,
                    str(warm), str(reps), "1" if inp == "sorted" else "0", str(bounds),
                ]
                # Timeout: the cap, plus a small floor so a hang is cut off.
                to = max(time_cap, 30.0)
                log(f"======= {name} Ckmeans k={k} input={inp} "
                    f"est_rss={est_rss:.1f}GiB est_fit={est_fit:.1f}s timeout={to:.0f}s =======")
                t0 = time.time()
                rc, out, err, peak, timed_out = run_wait4(cmd, timeout=to, env=env)
                elapsed = time.time() - t0
                with log_path.open("a") as f:
                    f.write(f"CMD {' '.join(cmd)} rc={rc} elapsed={elapsed:.3f}\n")
                    f.write(out)
                    if err:
                        f.write(err if err.endswith("\n") else err + "\n")
                kv = parse_cell(out)
                if rc == 124 or timed_out or not kv:
                    status = "timeout" if (rc == 124 or timed_out) else "crash"
                    row = dict(
                        dataset=name, method="Ckmeans.1d.dp", ckmeans_method=ckmeans_method,
                        k=k, n=n, D=D, input=inp, status=status,
                        sort_s=sort_s, exact_sse=exact, exact_sse_source=exact_src,
                        warm=warm, reps=reps,
                        peak_rss_mb="" if peak is None else f"{peak:.1f}",
                        est_rss_gib=f"{est_rss:.3f}", est_fit_s=f"{est_fit:.3f}",
                        elapsed_s=f"{elapsed:.3f}",
                        note=(err or out)[-500:],
                    )
                    append_row(csv_path, row)
                    rows.append(row)
                    log(f"{status.upper()} {name} k={k} input={inp} after {elapsed:.1f}s")
                    continue
                if kv.get("status") == "skip_k_ge_D":
                    row = {k2: kv.get(k2, "") for k2 in CSV_FIELDS}
                    row.update(dataset=name, method="Ckmeans.1d.dp",
                               ckmeans_method=kv.get("ckmeans_method", ckmeans_method),
                               k=k, n=kv.get("n", n), D=kv.get("D", D), input=inp,
                               status="skip_k_ge_D", sort_s=sort_s,
                               elapsed_s=f"{elapsed:.3f}")
                    append_row(csv_path, row)
                    rows.append(row)
                    log(f"skip_k_ge_D {name} k={k}")
                    continue

                sse = ""
                if bounds.is_file():
                    try:
                        sse = score_bounds(score, sorted_path, bounds)
                    except Exception as e:
                        log(f"warning: scoring failed for {name} k={k} input={inp}: {e}")

                ree = ""
                if sse and exact:
                    try:
                        ex = float(exact)
                        if ex > 0:
                            ree = f"{(float(sse) - ex) / ex * 100.0:.6g}"
                    except ValueError:
                        pass

                row = {k2: kv.get(k2, "") for k2 in CSV_FIELDS}
                row.update(
                    dataset=name,
                    method="Ckmeans.1d.dp",
                    ckmeans_method=kv.get("ckmeans_method", ckmeans_method),
                    version=kv.get("version", ""),
                    k=k,
                    n=kv.get("n", n),
                    D=kv.get("D", D),
                    input=inp,
                    status=kv.get("status", "ok"),
                    sort_s=sort_s,
                    post_s=kv.get("post_s", ""),
                    time_min_s=kv.get("time_min_s", ""),
                    time_max_s=kv.get("time_max_s", ""),
                    sse=sse,
                    tot_withinss=kv.get("tot_withinss", ""),
                    exact_sse=exact,
                    ree_pct=ree,
                    exact_sse_source=exact_src,
                    k_eff=kv.get("k_eff", ""),
                    warm=kv.get("warm", warm),
                    reps=kv.get("reps", reps),
                    peak_rss_mb="" if peak is None else f"{peak:.1f}",
                    est_rss_gib=f"{est_rss:.3f}",
                    est_fit_s=f"{est_fit:.3f}",
                    elapsed_s=f"{elapsed:.3f}",
                    note="",
                )
                append_row(csv_path, row)
                rows.append(row)
                log(f"CELL dataset={name} method=Ckmeans.1d.dp k={k} input={inp} "
                    f"status={row['status']} post_s={row['post_s']} sse={sse} "
                    f"ree_pct={ree} peak_rss_mb={row['peak_rss_mb']}")

    rows = load_csv(csv_path)
    write_summary(rows, OUT / "summary.md", smoke, mem_cap, time_cap, ckmeans_method)
    zip_path = PKG / "out" / "ckmeans_matched_mac_out.zip"
    if zip_path.exists():
        zip_path.unlink()
    subprocess.run(
        ["zip", "-r", "out/ckmeans_matched_mac_out.zip", "out/ckmeans_matched", "out/env.txt"],
        cwd=str(PKG), check=True,
    )
    log(f"Done in {time.time()-t_all:.1f}s. Results: out/ckmeans_matched/summary.md  "
        f"and  out/ckmeans_matched_mac_out.zip")


if __name__ == "__main__":
    main()
