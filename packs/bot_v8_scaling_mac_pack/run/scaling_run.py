#!/usr/bin/env python3
"""Orchestrate the v8 scaling run. Invoked by run/scaling.sh; not run by hand."""
from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
DATA = PKG / "data"
WORK = PKG / "work"
OUT = PKG / "out" / "scaling"
BUILD = PKG / "build"

KS = (64, 256, 512)
METHODS = ("wilber", "tuned", "sig5")  # wilber first so SSE* is known
FAMILIES = ("uniform", "normal", "lognormal")

AUDIO_GZ_NAME = "dharini-abhogi.txt.gz"
AUDIO_GZ_SHA = "3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27"
AUDIO_TXT_SHA = "b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e"
AUDIO_N = 8_168_256
T1_FILE = "gpt2-large_h17_attn_c_attn_weight.f32le"
T1_SHA = "04c28282e96a10b1375643345351f09dbf3874aaf22a63d87d968fa088b42a72"
T1_N = 4_915_200
T2_FILE = "gpt2-large_h35_mlp_c_proj_weight.f32le"
T2_SHA = "283c5aefc44e1c36ec8862e47589dec3fda4465f7c43618aa6dedd614f503836"
T2_N = 6_553_600

CSV_FIELDS = [
    "dataset", "method", "k", "M", "n", "D", "Dprime", "status",
    "sort_s", "post_s", "dp_s", "group_s", "atoms_s", "collapse_s", "refine_s",
    "sse", "exact_sse", "ree_pct", "coarse_sse", "coarse_ree_pct", "lib_cost",
    "k_eff", "atoms", "G", "collapse_passes", "collapse_merges", "merges_per_pass",
    "refine_sweeps", "refine_moves", "reps", "peak_rss_mb", "elapsed_s", "note",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_cell(line: str) -> dict | None:
    if not line.startswith("CELL "):
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


def done(rows: list[dict], dataset: str, method: str, k: int) -> bool:
    for r in rows:
        if r.get("dataset") == dataset and r.get("method") == method and int(r.get("k") or 0) == k:
            if r.get("status") in ("ok", "skip_k_ge_D"):
                return True
    return False


def default_M(k: int) -> int:
    return max(4096, 32 * k)


def wilber_timeout(n: int) -> float:
    # Box: ~15 s compute at N=8e6, any k. 8x margin, floor 45 s so a hang is cut off.
    return max(45.0, 120.0 * (n / 8_000_000))


def find_bin(name: str) -> Path:
    for p in (BUILD / name, BUILD / "Release" / name):
        if p.is_file() and os.access(p, os.X_OK):
            return p
    sys.exit(f"error: missing binary {name} under {BUILD}")


def run_cell(cmd: list[str], timeout: float | None) -> tuple[int, str, str]:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        return 124, out, err + f"\ntimeout after {timeout:.0f}s"


def parse_manifest(path: Path) -> dict[str, str]:
    rec = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        # family N sha256 file
        if len(parts) >= 4:
            rec[parts[3]] = parts[2]
        elif len(parts) == 2:
            rec[parts[1] if "/" in parts[1] or parts[1].endswith(".f64") else parts[0]] = parts[0]
    return rec


def generate_synth(python: str, exp: str) -> None:
    out = DATA / "synth"
    out.mkdir(parents=True, exist_ok=True)
    expect = parse_manifest(out / "MANIFEST.sha256")
    if not expect:
        sys.exit("error: missing data/synth/MANIFEST.sha256 in the pack")
    # generate only missing or mismatched
    need = []
    for fname, digest in expect.items():
        if not fname.endswith(".f64"):
            continue
        # honour --exp
        m = re.search(r"_n(\d+)\.f64$", fname)
        if m and exp:
            n = int(m.group(1))
            allowed = {1 << e for e in parse_exp(exp)}
            if n not in allowed:
                continue
        path = out / fname
        if path.is_file() and sha256_file(path) == digest:
            continue
        need.append(fname)
    # Keep the shipped full-grid manifest even if this run only generates a
    # subset (smoke uses 16-18). gen_synth.py would otherwise overwrite it.
    man_path = out / "MANIFEST.sha256"
    json_path = out / "MANIFEST.json"
    man_bytes = man_path.read_bytes()
    json_bytes = json_path.read_bytes() if json_path.is_file() else None
    if need:
        log(f"generating synthetic data (--exp {exp}) with {python} ...")
        r = subprocess.run([python, str(PKG / "run" / "gen_synth.py"),
                            "--exp", exp, "--out", str(out),
                            "--no-write-manifest"], check=False)
        man_path.write_bytes(man_bytes)
        if json_bytes is not None:
            json_path.write_bytes(json_bytes)
        if r.returncode != 0:
            sys.exit("error: gen_synth.py failed (need NumPy 1.17+, PCG64/default_rng)")
    else:
        man_path.write_bytes(man_bytes)
        if json_bytes is not None:
            json_path.write_bytes(json_bytes)
    # verify (only the requested exponents)
    allowed = {1 << e for e in parse_exp(exp)}
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
        # exp/log in normal/lognormal can differ in the last bits across
        # platforms (macOS libm vs Linux glibc). Seeds are the same; record the
        # local hashes and continue.
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


def parse_exp(s: str) -> list[int]:
    s = s.strip()
    if "-" in s and "," not in s:
        a, b = s.split("-", 1)
        lo, hi = int(a), int(b)
        return list(range(lo, hi + 1))
    return [int(x) for x in s.split(",") if x.strip()]


def prepare_optional(python: str, scale: Path, f32: Path, txt64: Path) -> list[tuple[str, Path]]:
    """Return (dataset, raw_f64_path) for audio and GPT-2 if present or fetchable."""
    found = []
    # audio
    gz = DATA / AUDIO_GZ_NAME
    txt = DATA / "dharini-abhogi.txt"
    work_txt = WORK / "dharini-abhogi.txt"
    audio_f64 = WORK / "audio.f64"
    if audio_f64.is_file() and audio_f64.stat().st_size == AUDIO_N * 8:
        found.append(("audio", audio_f64))
    elif gz.is_file() or txt.is_file():
        if gz.is_file():
            subprocess.run(["gunzip", "-c", str(gz)], check=True, stdout=work_txt.open("wb"))
        else:
            shutil.copy(txt, work_txt)
        got = sha256_file(work_txt)
        if got != AUDIO_TXT_SHA:
            log(f"warning: audio text SHA-256 {got} != {AUDIO_TXT_SHA}; skipping audio")
        else:
            subprocess.run([str(txt64), str(work_txt), str(audio_f64)], check=True)
            if audio_f64.stat().st_size != AUDIO_N * 8:
                sys.exit("error: audio f64 has the wrong size")
            found.append(("audio", audio_f64))
            log("audio ready")
    else:
        log("audio not in data/; skipping (place dharini-abhogi.txt.gz as usual to include it)")

    # GPT-2: fetch if network works, or accept copies. Skip if neither.
    have = all((DATA / f).is_file() for f in (T1_FILE, T2_FILE))
    if not have:
        log("GPT-2 tensors not in data/; trying range fetch ...")
        r = subprocess.run([python, str(PKG / "run" / "gpt2_fetch.py")], check=False)
        have = r.returncode == 0 and all((DATA / f).is_file() for f in (T1_FILE, T2_FILE))
        if not have:
            log("GPT-2 not available; skipping h17/h35 (copy the two .f32le files into data/ or allow network)")
    if have:
        for fname, digest, n, name in (
            (T1_FILE, T1_SHA, T1_N, "h17_c_attn"),
            (T2_FILE, T2_SHA, T2_N, "h35_c_proj"),
        ):
            path = DATA / fname
            got = sha256_file(path)
            if got != digest:
                log(f"warning: {fname} SHA-256 {got} != {digest}; skipping {name}")
                continue
            outp = WORK / f"{name}.f64"
            if not (outp.is_file() and outp.stat().st_size == n * 8):
                subprocess.run([str(f32), str(path), str(outp)], check=True)
            found.append((name, outp))
        if any(n in {x[0] for x in found} for n in ("h17_c_attn", "h35_c_proj")):
            log("GPT-2 tensors ready")
    return found


def write_summary(rows: list[dict], path: Path, smoke: bool) -> None:
    # reload
    by = {}
    sort_s = {}
    for r in rows:
        if r.get("method") == "sort":
            sort_s[r["dataset"]] = r
            continue
        by[(r["dataset"], r.get("method"), int(r.get("k") or 0))] = r

    def fnum(r, k, nd=3):
        try:
            return f"{float(r.get(k) or 0):.{nd}f}"
        except ValueError:
            return r.get(k) or ""

    lines = []
    lines.append("# Scaling: lex1dx tuned / sig5 (always 5-sig) vs Grønlund wilber-interp")
    lines.append("")
    if smoke:
        lines.append("SMOKE: synthetic N = 2^16..2^18, all three families, k ∈ {64,256,512}, "
                     "one warm-up and one timed rep. Times are that single rep, not a median.")
    else:
        lines.append("One warm-up and 3 timed reps; times are medians, in seconds. "
                     "**Comparison is post-sort only**: every method is given the same array, "
                     "already sorted with `std::sort` (that sort is timed separately and logged "
                     "as `sort_s`; it is not added to `post_s`). "
                     "lex1dx **tuned** uses full-precision atoms, tuned Collapse, weighted wilber DP, Refine. "
                     "lex1dx **sig5** always groups to 5 significant digits (no D′ guard), then the same pipeline; D and D′ are logged so a rule can be derived later. "
                     "M = max(4096, 32k). Exact SSE* is the two-pass long-double SSE of unweighted Grønlund wilber-interp on the sorted points. "
                     "REE is (SSE − SSE*) / SSE* × 100. Peak RSS is the cell process (`getrusage`).")
    lines.append("")
    lines.append("If a wilber-interp cell timed out or crashed, its row has `status=timeout` or `status=crash` and lex REE is left blank.")
    lines.append("")
    # sort reference
    if sort_s:
        lines.append("## Sort reference (`std::sort`, not in post-sort)")
        lines.append("")
        lines.append("| dataset | n | D | sort s | peak RSS MiB |")
        lines.append("|---|---:|---:|---:|---:|")
        for ds, r in sorted(sort_s.items(), key=lambda kv: (kv[0])):
            lines.append(f"| {ds} | {r.get('n','')} | {r.get('D','')} | {fnum(r,'sort_s',3)} | {fnum(r,'peak_rss_mb',1)} |")
        lines.append("")

    lines.append("## Post-sort")
    lines.append("")
    lines.append("| dataset | k | M | method | status | D | D' | group s | collapse s | DP s | refine s | post-sort s | vs wilber | SSE | REE % | passes | merges | refine moves | RSS MiB |")
    lines.append("|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    datasets = sorted({k[0] for k in by}, key=lambda d: (0, int(d.split("_n")[-1])) if "_n" in d else (1, d))
    for ds in datasets:
        for k in KS:
            w = by.get((ds, "wilber", k))
            wpost = float(w["post_s"]) if w and w.get("status") == "ok" and w.get("post_s") else None
            for meth in METHODS:
                r = by.get((ds, meth, k))
                if not r:
                    continue
                post = float(r["post_s"]) if r.get("post_s") else float("nan")
                sp = f"{wpost / post:.2f}×" if wpost and post and post > 0 else ""
                ree = r.get("ree_pct") or ""
                try:
                    ree = f"{float(ree):.5g}" if ree not in ("", "nan", "NAN") else ""
                except ValueError:
                    pass
                lines.append(
                    f"| {ds} | {k} | {r.get('M') or default_M(k)} | {meth} | {r.get('status')} | "
                    f"{r.get('D')} | {r.get('Dprime')} | {fnum(r,'group_s')} | {fnum(r,'collapse_s')} | "
                    f"{fnum(r,'dp_s')} | {fnum(r,'refine_s')} | {fnum(r,'post_s')} | {sp} | "
                    f"{r.get('sse') or ''} | {ree} | {r.get('collapse_passes') or ''} | "
                    f"{r.get('collapse_merges') or ''} | {r.get('refine_moves') or ''} | "
                    f"{fnum(r,'peak_rss_mb',1)} |"
                )
    lines.append("")
    lines.append("See `out/scaling/runs.csv` and `out/scaling/run.log`. Re-run the same command to resume; finished (status=ok) cells are skipped.")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    smoke = "--smoke" in sys.argv
    python = os.environ.get("PYTHON_BIN") or sys.executable
    cell = find_bin("scale_cell")
    f32 = find_bin("f32_to_f64")
    txt64 = find_bin("txt_to_f64")
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "runs.csv"
    log_path = OUT / "run.log"
    rows = load_csv(csv_path)

    exp = "16-18" if smoke else "16-23"
    reps = 1 if smoke else 3
    generate_synth(python, exp)

    datasets: list[tuple[str, Path]] = []
    allowed = {1 << e for e in parse_exp(exp)}
    for fam in FAMILIES:
        for n in sorted(allowed):
            p = DATA / "synth" / f"{fam}_n{n}.f64"
            if p.is_file():
                datasets.append((f"{fam}_n{n}", p))
    if not smoke:
        datasets.extend(prepare_optional(python, cell, f32, txt64))
    if not datasets:
        sys.exit("error: no datasets to run")

    log(f"datasets={len(datasets)} ks={list(KS)} reps={reps} smoke={int(smoke)}")
    t_all = time.time()

    for name, raw in datasets:
        n = raw.stat().st_size // 8
        sorted_path = WORK / f"{name}.sorted.f64"
        if not done(rows, name, "sort", 0) or not sorted_path.is_file():
            log(f"======= sort {name} n={n} =======")
            rc, out, err = run_cell(
                [str(cell), "sort", str(raw), str(sorted_path), name, str(reps)],
                timeout=max(30.0, 60.0 * (n / 8_000_000)),
            )
            (log_path.open("a")).write(out + err + f"\nrc={rc}\n")
            cell_line = next((ln for ln in out.splitlines() if ln.startswith("CELL ")), None)
            if rc != 0 or not cell_line:
                row = dict(dataset=name, method="sort", k=0, n=n, status="timeout" if rc == 124 else "crash",
                           note=(err or out)[-400:], reps=reps)
                append_row(csv_path, row)
                rows.append(row)
                log(f"SORT FAILED {name} rc={rc}")
                continue
            kv = parse_cell(cell_line)
            row = {k: kv.get(k, "") for k in CSV_FIELDS}
            row["elapsed_s"] = ""
            append_row(csv_path, row)
            rows.append(row)
            log(cell_line)
        else:
            log(f"skip sort {name} (already done)")

        if not sorted_path.is_file():
            continue

        exact = {}  # k -> sse
        for r in rows:
            if r.get("dataset") == name and r.get("method") == "wilber" and r.get("status") == "ok":
                try:
                    exact[int(r["k"])] = r["sse"]
                except (KeyError, ValueError):
                    pass

        for k in KS:
            for meth in METHODS:
                if done(rows, name, meth, k):
                    log(f"skip {name} {meth} k={k}")
                    continue
                ex = exact.get(k, "")
                cmd = [str(cell), meth, str(sorted_path), name, str(k), str(reps)]
                if meth != "wilber" and ex:
                    cmd.append(str(ex))
                to = wilber_timeout(n) if meth == "wilber" else max(60.0, 30.0 * (n / 8_000_000) + k / 64.0)
                log(f"======= {name} {meth} k={k} timeout={to:.0f}s =======")
                t0 = time.time()
                rc, out, err = run_cell(cmd, timeout=to)
                elapsed = time.time() - t0
                with log_path.open("a") as f:
                    f.write(f"CMD {' '.join(cmd)} rc={rc} elapsed={elapsed:.3f}\n")
                    f.write(out)
                    if err:
                        f.write(err)
                        if not err.endswith("\n"):
                            f.write("\n")
                cell_line = next((ln for ln in out.splitlines() if ln.startswith("CELL ")), None)
                if rc == 124 or not cell_line:
                    status = "timeout" if rc == 124 else "crash"
                    row = dict(dataset=name, method=meth, k=k, M=default_M(k) if meth != "wilber" else 0,
                               n=n, status=status, elapsed_s=f"{elapsed:.3f}",
                               note=(err or out)[-500:], reps=reps)
                    append_row(csv_path, row)
                    rows.append(row)
                    log(f"{status.upper()} {name} {meth} k={k} after {elapsed:.1f}s")
                    continue
                kv = parse_cell(cell_line)
                # attach sort_s from the sort cell
                sr = sort_s_row = next((r for r in rows if r.get("dataset") == name and r.get("method") == "sort"), None)
                if sr and (not kv.get("sort_s") or kv.get("sort_s") in ("0", "0.000000")):
                    kv["sort_s"] = sr.get("sort_s", "")
                if meth != "wilber" and ex and not kv.get("exact_sse"):
                    kv["exact_sse"] = ex
                if meth == "wilber" and kv.get("status") == "ok" and kv.get("sse"):
                    exact[k] = kv["sse"]
                row = {k2: kv.get(k2, "") for k2 in CSV_FIELDS}
                row["elapsed_s"] = f"{elapsed:.3f}"
                append_row(csv_path, row)
                rows.append(row)
                log(cell_line)

    rows = load_csv(csv_path)
    write_summary(rows, OUT / "summary.md", smoke)
    zip_path = PKG / "out" / "scaling_mac_out.zip"
    if zip_path.exists():
        zip_path.unlink()
    subprocess.run(
        ["zip", "-r", "out/scaling_mac_out.zip", "out", "-x", "out/scaling_mac_out.zip"],
        cwd=str(PKG), check=True,
    )
    log(f"Done in {time.time()-t_all:.1f}s. Results: out/scaling/summary.md  and  out/scaling_mac_out.zip")


if __name__ == "__main__":
    main()
