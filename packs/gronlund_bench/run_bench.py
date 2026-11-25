#!/usr/bin/env python3
"""Run correctness and small-scale benches. Writes CSV + per-cell logs."""
from __future__ import annotations

import csv
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/workspace/lexicon/work/gronlund_bench")
DRIVER = ROOT / "gronlund_driver"
LEX = ROOT / "lex_build" / "lex_cli"
DATA = ROOT / "data"
LOGS = ROOT / "logs"
OUT = ROOT / "results"
LOGS.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

# overall scale budget (seconds), leaving a little room to write the report
SCALE_BUDGET_S = 18 * 60
CORRECT_TIMEOUT = 60
SCALE_TIMEOUT = 180  # per cell; skip if remaining budget is below this for large cells

GRONLUND_ALGS = [
    ("linear", "report"),
    ("monotone", "report"),
    ("slow", "report"),
    ("hirc", "report"),
    ("wilber-binary", "report"),
    ("wilber-interp", "report"),
]
# scale: cost-only and report for the two-row DP methods; no slow
SCALE_ALGS = [
    ("linear", "compute"),
    ("linear", "report"),
    ("monotone", "compute"),
    ("monotone", "report"),
    ("hirc", "report"),
    ("wilber-binary", "report"),
    ("wilber-interp", "report"),
]

RESULT_RE = re.compile(r"^RESULT\s+(.*)$", re.M)
SSE_RE = re.compile(r"^sse=([0-9.eE+-]+)$", re.M)


def parse_kv(s: str) -> dict:
    out = {}
    for part in s.split():
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def parse_time_file(path: Path) -> dict:
    text = path.read_text() if path.exists() else ""
    d = {"time_text": text}
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([0-9.:]+)", text)
    if m:
        t = m.group(1)
        if t.count(":") == 2:
            h, mn, s = t.split(":")
            d["elapsed_s"] = int(h) * 3600 + int(mn) * 60 + float(s)
        else:
            mn, s = t.split(":")
            d["elapsed_s"] = int(mn) * 60 + float(s)
    m = re.search(r"Maximum resident set size \(kbytes\): ([0-9]+)", text)
    if m:
        d["max_rss_kb"] = int(m.group(1))
        d["max_rss_mb"] = int(m.group(1)) / 1024.0
    m = re.search(r"Exit status: ([0-9]+)", text)
    if m:
        d["time_exit"] = int(m.group(1))
    m = re.search(r"Command terminated by signal (\d+)", text)
    if m:
        d["signal"] = int(m.group(1))
    return d


def run_cmd(cmd, log_prefix: Path, timeout: int) -> dict:
    stdout_p = Path(str(log_prefix) + ".stdout")
    stderr_p = Path(str(log_prefix) + ".stderr")
    time_p = Path(str(log_prefix) + ".time")
    wrapper = [
        "/usr/bin/time", "-v", "-o", str(time_p),
        "timeout", f"--signal=KILL", str(timeout),
        *cmd,
    ]
    t0 = time.time()
    try:
        p = subprocess.run(
            wrapper,
            stdout=open(stdout_p, "w"),
            stderr=open(stderr_p, "w"),
            cwd=str(ROOT),
        )
        rc = p.returncode
    except Exception as e:
        rc = -1
        stderr_p.write_text(str(e))
    wall = time.time() - t0
    stdout = stdout_p.read_text() if stdout_p.exists() else ""
    stderr = stderr_p.read_text() if stderr_p.exists() else ""
    td = parse_time_file(time_p)
    rec = {
        "rc": rc,
        "wall_s": wall,
        "stdout": stdout,
        "stderr": stderr,
        **{k: v for k, v in td.items() if k != "time_text"},
    }
    # timeout(1) returns 124; with KILL, 137
    if rc in (124, 137):
        rec["status"] = "timeout"
    elif rc != 0:
        sig = rec.get("signal")
        rec["status"] = f"crash rc={rc}" + (f" signal={sig}" if sig else "")
        # common: 134 SIGABRT (assert), 139 SIGSEGV, 136 SIGFPE
        if rc == 134 or sig == 6:
            rec["status"] = "abort (SIGABRT, likely assert)"
        elif rc == 139 or sig == 11:
            rec["status"] = "segfault (SIGSEGV)"
        elif rc == 136 or sig == 8:
            rec["status"] = "SIGFPE"
    else:
        rec["status"] = "ok"
    return rec


def run_gronlund(dataset: Path, k: int, alg: str, mode: str, timeout: int, tag: str) -> dict:
    log = LOGS / f"{tag}_{dataset.stem}_k{k}_{alg}_{mode}"
    cmd = [str(DRIVER), "--alg", alg, "--k", str(k), "--mode", mode, "--out", str(log) + ".result", str(dataset)]
    rec = run_cmd(cmd, log, timeout)
    rec.update(dataset=dataset.name, n=None, k=k, method=f"{alg}:{mode}", family="gronlund")
    m = RESULT_RE.search(rec["stdout"])
    if m:
        kv = parse_kv(m.group(1))
        rec.update(kv)
        rec["sse"] = float(kv["sse"]) if "sse" in kv else None
        rec["n"] = int(kv["n"]) if "n" in kv else None
        rec["driver_s"] = float(kv["time_s"]) if "time_s" in kv else None
        rec["n_centers"] = int(kv["n_centers"]) if "n_centers" in kv else None
    return rec


def run_lex(dataset: Path, k: int, preset: str, timeout: int, tag: str, coarse=None) -> dict:
    log = LOGS / (f"{tag}_{dataset.stem}_k{k}_lex_{preset}" + (f"_M{coarse}" if coarse is not None else ""))
    cmd = [str(LEX), str(dataset), str(k), "--preset", preset]
    if coarse is not None:
        cmd += ["--coarse", str(coarse)]
    rec = run_cmd(cmd, log, timeout)
    rec.update(dataset=dataset.name, n=None, k=k, method=f"lex_{preset}" + (f"_M{coarse}" if coarse is not None else ""), family="lex")
    m = SSE_RE.search(rec["stdout"])
    if m:
        rec["sse"] = float(m.group(1))
    # n from wc if needed later
    return rec


def rel_diff(sse, sse_star):
    if sse is None or sse_star is None:
        return None
    if sse_star == 0.0:
        return 0.0 if sse == 0.0 else float("inf")
    return (sse - sse_star) / sse_star


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    t_all = time.time()
    # ---------- correctness ----------
    kinds = ["uniform", "normal", "lognormal", "ties"]
    ns = [1000, 10000]
    ks = [4, 16, 64]
    correct_rows = []
    for n in ns:
        for kind in kinds:
            ds = DATA / f"{kind}_n{n}.txt"
            for k in ks:
                ref = run_lex(ds, k, "exact", CORRECT_TIMEOUT, "corr")
                ref_row = {
                    "phase": "correctness",
                    "dataset": ds.name,
                    "kind": kind,
                    "n": n,
                    "k": k,
                    "method": ref["method"],
                    "status": ref["status"],
                    "sse": ref.get("sse"),
                    "sse_rel_vs_exact": 0.0 if ref.get("sse") is not None else None,
                    "n_centers": None,
                    "driver_s": None,
                    "elapsed_s": ref.get("elapsed_s"),
                    "max_rss_mb": ref.get("max_rss_mb"),
                    "rc": ref["rc"],
                    "error": (ref["stderr"] or "")[:400] if ref["status"] != "ok" else "",
                }
                correct_rows.append(ref_row)
                sse_star = ref.get("sse")
                algs = list(GRONLUND_ALGS)
                # slow is O(k n^2); keep it for all correctness cells but timeout will catch overruns
                for alg, mode in algs:
                    rec = run_gronlund(ds, k, alg, mode, CORRECT_TIMEOUT, "corr")
                    row = {
                        "phase": "correctness",
                        "dataset": ds.name,
                        "kind": kind,
                        "n": n,
                        "k": k,
                        "method": rec["method"],
                        "status": rec["status"],
                        "sse": rec.get("sse"),
                        "sse_rel_vs_exact": rel_diff(rec.get("sse"), sse_star),
                        "n_centers": rec.get("n_centers"),
                        "driver_s": rec.get("driver_s"),
                        "elapsed_s": rec.get("elapsed_s"),
                        "max_rss_mb": rec.get("max_rss_mb"),
                        "rc": rec["rc"],
                        "error": (rec["stderr"] or rec["stdout"][:400] if rec["status"] != "ok" else ""),
                    }
                    correct_rows.append(row)
                    print(f"CORR {ds.name} k={k} {alg}:{mode} status={row['status']} rel={row['sse_rel_vs_exact']}", flush=True)

    fields = ["phase","dataset","kind","n","k","method","status","sse","sse_rel_vs_exact","n_centers","driver_s","elapsed_s","max_rss_mb","rc","error"]
    write_csv(OUT / "correctness.csv", correct_rows, fields)

    # max |rel| among ok rows with finite rel, excluding the exact reference itself
    rels = []
    for r in correct_rows:
        if r["status"] == "ok" and r["method"] != "lex_exact" and r["sse_rel_vs_exact"] is not None:
            v = r["sse_rel_vs_exact"]
            if v == float("inf"):
                rels.append(("inf", r))
            else:
                rels.append((abs(v), r))
    finite = [a for a, _ in rels if a != "inf"]
    infs = [r for a, r in rels if a == "inf"]
    max_rel = max(finite) if finite else None
    (OUT / "correctness_summary.txt").write_text(
        f"n_rows={len(correct_rows)}\n"
        f"n_ok_compared={len(finite)+len(infs)}\n"
        f"n_infinite_rel (exact SSE=0 but method SSE!=0)={len(infs)}\n"
        f"max_abs_rel_sse={max_rel}\n"
        + ("\n".join(
            f"  inf: {r['dataset']} k={r['k']} {r['method']} sse={r['sse']} exact=0"
            for r in infs
        ) + "\n" if infs else "")
    )
    print("=== correctness done, max_abs_rel=", max_rel, "infs", len(infs), flush=True)

    # ---------- scale ----------
    scale_rows = []
    scale_start = time.time()
    scale_jobs = []
    for ds, n_hint, kind in [
        (DATA / "normal_n100000.txt", 100000, "normal"),
        (DATA / "normal_n1000000.txt", 1000000, "normal"),
        (DATA / "audio_prefix_n100000.txt", 100000, "audio_prefix"),
        (DATA / "audio_prefix_n1000000.txt", 1000000, "audio_prefix"),
    ]:
        for k in (64, 512):
            # LEX first (needed? not as reference at this scale; still time it)
            scale_jobs.append(("lex", ds, n_hint, kind, k, None, None))
            for alg, mode in SCALE_ALGS:
                scale_jobs.append(("gronlund", ds, n_hint, kind, k, alg, mode))

    skipped = 0
    for job in scale_jobs:
        used = time.time() - scale_start
        remaining = SCALE_BUDGET_S - used
        typ, ds, n_hint, kind, k, alg, mode = job
        # skip a cell if remaining budget is tight relative to problem size
        # (do not start N=1e6 k=512 if under 30s left, etc.)
        need = 5
        if n_hint >= 1000000 and k >= 512:
            need = 60
        elif n_hint >= 1000000:
            need = 30
        elif k >= 512:
            need = 20
        label = f"{ds.name} k={k} " + (f"{alg}:{mode}" if typ == "gronlund" else "lex_hybrid_M4096")
        if remaining < need:
            print(f"SKIP {label} remaining={remaining:.1f}s need={need}s", flush=True)
            scale_rows.append({
                "phase": "scale",
                "dataset": ds.name,
                "kind": kind,
                "n": n_hint,
                "k": k,
                "method": f"{alg}:{mode}" if typ == "gronlund" else "lex_hybrid_M4096",
                "status": f"skipped (remaining {remaining:.0f}s < {need}s budget headroom)",
                "sse": None,
                "n_centers": None,
                "driver_s": None,
                "elapsed_s": None,
                "max_rss_mb": None,
                "rc": None,
                "error": "",
            })
            skipped += 1
            continue
        if typ == "lex":
            rec = run_lex(ds, k, "hybrid", SCALE_TIMEOUT, "scale", coarse=4096)
            method = rec["method"]
        else:
            rec = run_gronlund(ds, k, alg, mode, SCALE_TIMEOUT, "scale")
            method = rec["method"]
        row = {
            "phase": "scale",
            "dataset": ds.name,
            "kind": kind,
            "n": rec.get("n") or n_hint,
            "k": k,
            "method": method,
            "status": rec["status"],
            "sse": rec.get("sse"),
            "n_centers": rec.get("n_centers"),
            "driver_s": rec.get("driver_s"),
            "elapsed_s": rec.get("elapsed_s"),
            "max_rss_mb": rec.get("max_rss_mb"),
            "rc": rec["rc"],
            "error": (rec.get("stderr") or "")[:400] if rec["status"] != "ok" else "",
        }
        scale_rows.append(row)
        print(f"SCALE {label} status={row['status']} elapsed={row['elapsed_s']} rss_mb={row['max_rss_mb']} sse={row['sse']}", flush=True)

    sfields = ["phase","dataset","kind","n","k","method","status","sse","n_centers","driver_s","elapsed_s","max_rss_mb","rc","error"]
    write_csv(OUT / "scale.csv", scale_rows, sfields)
    print(f"=== all done in {time.time()-t_all:.1f}s  scale_used={time.time()-scale_start:.1f}s skipped={skipped}", flush=True)


if __name__ == "__main__":
    main()
