#!/usr/bin/env python3
"""Orchestrate the v9 controls run. Invoked by run/controls.sh; not run by hand.

Reuses data handling from run/scaling_run.py (v8). Three experiments:
  A  baseline A/B: pristine vendor/kmeans1d wilber vs patched bench/weighted_wilber (unit weights)
  B  sig5exact control: 5-sig grouping + exact weighted wilber on all D' groups, no Collapse/Refine
     (sig5 default is re-run alongside for a same-session comparison)
  C  Collapse quantile 0.5 / 0.75 / 0.9 for tuned and sig5, normal 8M + audio, k=512
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scaling_run as S  # noqa: E402  (v8 helpers: synth generation, csv, cell parsing)

PKG = S.PKG
DATA, WORK, BUILD = S.DATA, S.WORK, S.BUILD
OUT = PKG / "out" / "controls"
REF_V8 = PKG / "ref" / "v8_runs.csv"

KS = (64, 256, 512)
QS = ("0.5", "0.75", "0.9")
FIELDS = S.CSV_FIELDS[:-2] + ["q", "experiment", "elapsed_s", "note"]
log = S.log


def done(rows, ds, meth, k, q):
    return any(r.get("dataset") == ds and r.get("method") == meth and int(r.get("k") or 0) == k
               and (r.get("q") or "") == q and r.get("status") in ("ok", "skip_k_ge_D") for r in rows)


def append(path: Path, row: dict) -> None:
    new = not path.is_file()
    with path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDS})


def prepare_audio(txt64: Path, smoke: bool) -> Path | None:
    """Same audio handling as v8 (data/dharini-abhogi.txt.gz or .txt, hash checked)."""
    audio_f64 = WORK / "audio.f64"
    if smoke:
        # Fake audio: integer-valued, tie-heavy, small. Never the real file.
        fake = WORK / "audio_fake.f64"
        if not fake.is_file():
            import numpy as np
            rng = np.random.Generator(np.random.PCG64(7))
            x = np.round(rng.normal(0, 3000, 65536)).astype("<f8")
            x.tofile(str(fake))
        return fake
    if audio_f64.is_file() and audio_f64.stat().st_size == S.AUDIO_N * 8:
        return audio_f64
    gz, txt = DATA / S.AUDIO_GZ_NAME, DATA / "dharini-abhogi.txt"
    if not (gz.is_file() or txt.is_file()):
        log("warning: audio not in data/; audio cells skipped (place dharini-abhogi.txt.gz in data/)")
        return None
    work_txt = WORK / "dharini-abhogi.txt"
    if gz.is_file():
        with work_txt.open("wb") as f:
            subprocess.run(["gunzip", "-c", str(gz)], check=True, stdout=f)
    else:
        shutil.copy(txt, work_txt)
    got = S.sha256_file(work_txt)
    if got != S.AUDIO_TXT_SHA:
        log(f"warning: audio text SHA-256 {got} != {S.AUDIO_TXT_SHA}; continuing anyway")
    subprocess.run([str(txt64), str(work_txt), str(audio_f64)], check=True)
    log("audio ready")
    return audio_f64


def cell_timeout(meth: str, n: int, k: int) -> float:
    if meth.startswith("wilber"):
        return max(60.0, 300.0 * (n / 8_000_000))  # v8 Mac: <= 12 s per run, 4 runs per cell
    return max(60.0, 60.0 * (n / 8_000_000) + k / 32.0)


def main() -> None:
    smoke = "--smoke" in sys.argv
    python = os.environ.get("PYTHON_BIN") or sys.executable
    cell, pristine = S.find_bin("scale_cell"), S.find_bin("pristine_cell")
    txt64 = S.find_bin("txt_to_f64")
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path, log_path = OUT / "runs.csv", OUT / "run.log"
    rows = S.load_csv(csv_path)

    exp = "16" if smoke else "23"
    reps = 1 if smoke else 3
    S.generate_synth(python, exp)  # hash mismatch -> warning only (v8 behaviour)
    n0 = 1 << int(exp)
    datasets = []
    for fam in S.FAMILIES:
        p = DATA / "synth" / f"{fam}_n{n0}.f64"
        if p.is_file():
            datasets.append((f"{fam}_n{n0}", p))
    a = prepare_audio(txt64, smoke)
    if a is not None:
        datasets.append(("audio_fake" if smoke else "audio", a))
    qsets = {f"normal_n{n0}", "audio_fake" if smoke else "audio"}

    # plan: (dataset, experiment, method, k, q)
    plan = []
    for name, _ in datasets:
        for k in KS:
            plan.append((name, "A", "wilber", k, ""))
            plan.append((name, "A", "wilber_pristine", k, ""))
            plan.append((name, "B", "sig5exact", k, ""))
            plan.append((name, "B", "sig5", k, ""))
        if name in qsets:
            for meth in ("tuned", "sig5"):
                for q in QS:
                    plan.append((name, "C", meth, 512, q))
    log(f"datasets={[d for d, _ in datasets]} cells={len(plan)} reps={reps} smoke={int(smoke)}")
    t_all = time.time()

    for name, raw in datasets:
        n = raw.stat().st_size // 8
        sp = WORK / f"{name}.sorted.f64"
        if not done(rows, name, "sort", 0, "") or not sp.is_file():
            rc, out, err = S.run_cell([str(cell), "sort", str(raw), str(sp), name, str(reps)],
                                      timeout=max(60.0, 120.0 * n / 8_000_000))
            with log_path.open("a") as f:
                f.write(out + err + f"\nrc={rc}\n")
            line = next((ln for ln in out.splitlines() if ln.startswith("CELL ")), None)
            if rc != 0 or not line:
                append(csv_path, dict(dataset=name, method="sort", k=0, n=n, status="crash", note=(err or out)[-400:]))
                log(f"SORT FAILED {name} rc={rc}")
                continue
            kv = S.parse_cell(line)
            kv["experiment"] = "sort"
            append(csv_path, kv)
            rows.append(kv)
            log(line)
        sort_s = next((r.get("sort_s", "") for r in rows if r.get("dataset") == name and r.get("method") == "sort"), "")
        exact = {int(r["k"]): r["sse"] for r in rows
                 if r.get("dataset") == name and r.get("method") == "wilber" and r.get("status") == "ok"}

        for ds, xp, meth, k, q in plan:
            if ds != name:
                continue
            if done(rows, name, meth, k, q):
                log(f"skip {name} {meth} k={k} q={q or '-'} (done)")
                continue
            binp = pristine if meth == "wilber_pristine" else cell
            cmd = [str(binp), "wilber" if meth == "wilber_pristine" else meth, str(sp), name, str(k), str(reps)]
            if not meth.startswith("wilber") and k in exact:
                cmd.append(str(exact[k]))
            if q:
                cmd.append(f"--q={q}")
            to = cell_timeout(meth, n, k)
            log(f"======= [{xp}] {name} {meth} k={k} q={q or '-'} =======")
            t0 = time.time()
            rc, out, err = S.run_cell(cmd, timeout=to)
            el = time.time() - t0
            with log_path.open("a") as f:
                f.write(f"CMD {' '.join(cmd)} rc={rc} elapsed={el:.3f}\n{out}{err}\n")
            line = next((ln for ln in out.splitlines() if ln.startswith("CELL ")), None)
            if rc == 124 or not line:
                st = "timeout" if rc == 124 else "crash"
                # Logged, not truncated: a timeout is a lower bound; rerun with a longer limit.
                row = dict(dataset=name, method=meth, k=k, q=q, n=n, status=st, experiment=xp,
                           elapsed_s=f"{el:.3f}", note=f"limit={to:.0f}s " + (err or out)[-400:], reps=reps)
                append(csv_path, row)
                rows.append(row)
                log(f"{st.upper()} {name} {meth} k={k} after {el:.1f}s (limit {to:.0f}s)")
                continue
            kv = S.parse_cell(line)
            kv["method"] = meth
            kv["q"] = q
            kv["experiment"] = xp
            kv["sort_s"] = sort_s
            kv["elapsed_s"] = f"{el:.3f}"
            if meth == "wilber" and kv.get("status") == "ok":
                exact[k] = kv["sse"]
            if k in exact:
                kv["exact_sse"] = exact[k]
                try:
                    kv["ree_pct"] = "%.6g" % ((float(kv["sse"]) - float(exact[k])) / float(exact[k]) * 100.0)
                except (KeyError, ValueError, ZeroDivisionError):
                    pass
            append(csv_path, kv)
            rows.append(kv)
            log(line)

    write_summary(S.load_csv(csv_path), OUT / "summary.md", smoke, time.time() - t_all)
    zp = PKG / "out" / "controls_mac_out.zip"
    if zp.exists():
        zp.unlink()
    subprocess.run(["zip", "-qr", "out/controls_mac_out.zip", "out", "-x", "out/controls_mac_out.zip"],
                   cwd=str(PKG), check=False)
    log(f"Done in {time.time()-t_all:.1f}s. Results: out/controls/summary.md, runs.csv, run.log, out/env.txt "
        f"-> send out/controls_mac_out.zip")


def f(r, key, nd=3):
    try:
        return f"{float(r[key]):.{nd}f}"
    except (KeyError, ValueError, TypeError):
        return "—"


def write_summary(rows, path: Path, smoke: bool, wall: float) -> None:
    v8 = {}
    if REF_V8.is_file():
        for r in S.load_csv(REF_V8):
            v8[(r["dataset"], r["method"], int(r.get("k") or 0))] = r
    by = {(r["dataset"], r["method"], int(r.get("k") or 0), r.get("q") or ""): r
          for r in rows if r.get("status") == "ok"}
    dsets = []
    for r in rows:
        if r["dataset"] not in dsets:
            dsets.append(r["dataset"])
    L = [f"# v9 controls{' (SMOKE)' if smoke else ''}", "",
         "Post-sort timing only (std::sort excluded, shown as sort_s). Times are medians of the timed reps. "
         "All SSE is on the raw sorted points; REE is vs Grønlund SSE* (patched wilber, unit weights).", ""]
    bad = [r for r in rows if r.get("status") not in ("ok", "skip_k_ge_D")]
    L.append(f"Cells: {len(rows)} rows, {len(bad)} not ok" + (": " + ", ".join(
        f"{r['dataset']} {r['method']} k={r['k']} {r['status']}" for r in bad) if bad else "") + ".")
    L.append("")

    L += ["## A. Baseline A/B: pristine vs patched Grønlund wilber-interp", "",
          "| dataset | k | sort_s | pristine post_s | patched post_s | patched/pristine | SSE identical |",
          "|---|---:|---:|---:|---:|---:|:---:|"]
    for d in dsets:
        for k in KS:
            a, b = by.get((d, "wilber_pristine", k, "")), by.get((d, "wilber", k, ""))
            if not a or not b:
                continue
            ratio = float(b["post_s"]) / float(a["post_s"]) if float(a["post_s"]) > 0 else float("nan")
            L.append(f"| {d} | {k} | {f(b,'sort_s')} | {f(a,'post_s')} | {f(b,'post_s')} | {ratio:.3f} | "
                     f"{'yes' if a['sse'] == b['sse'] else 'NO ' + a['sse'] + ' vs ' + b['sse']} |")
    L.append("")

    L += ["## B. sig5exact control (5-sig rounding + exact weighted DP on all D′ groups; no Collapse, no Refine)", "",
          "| dataset | k | D′ | sig5exact post_s | sig5exact REE % | sig5exact RSS MB | sig5 post_s | sig5 REE % | "
          "sig5 RSS MB | v8 sig5 post_s | wilber post_s |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for d in dsets:
        for k in KS:
            x, s, w = by.get((d, "sig5exact", k, "")), by.get((d, "sig5", k, "")), by.get((d, "wilber", k, ""))
            if not x:
                continue
            s = s or {}
            w = w or {}
            v = v8.get((d, "sig5", k), {})
            L.append(f"| {d} | {k} | {x.get('Dprime','')} | {f(x,'post_s')} | {f(x,'ree_pct',5)} | {f(x,'peak_rss_mb',1)} | "
                     f"{f(s,'post_s')} | {f(s,'ree_pct',5)} | {f(s,'peak_rss_mb',1)} | {f(v,'post_s')} | {f(w,'post_s')} |")
    L.append("")

    L += ["## C. Collapse quantile sensitivity (k = 512; 0.75 is the default)", "",
          "| dataset | method | q | post_s | collapse_s | REE % | collapse passes | G |",
          "|---|---|---:|---:|---:|---:|---:|---:|"]
    for d in dsets:
        for m in ("tuned", "sig5"):
            for q in QS:
                r = by.get((d, m, 512, q))
                if r:
                    L.append(f"| {d} | {m} | {q} | {f(r,'post_s')} | {f(r,'collapse_s')} | {f(r,'ree_pct',5)} | "
                             f"{r.get('collapse_passes','')} | {r.get('G','')} |")
    L.append("")

    L += ["## Reproducibility vs v8 (same dataset, method, k, default q)", ""]
    L.append("tuned and wilber should match v8 (same pipeline). sig5 is expected to differ: "
             "v10 Refines on raw distinct values and raises grouping precision when D' < k.")
    chk = []
    for (d, m, k, q), r in by.items():
        mm = m if m in ("wilber", "tuned") else None
        if mm is None or q not in ("", "0.75"):
            continue
        ref = v8.get((d, mm, k))
        if ref and ref.get("sse"):
            chk.append((d, m, k, q, r["sse"] == ref["sse"], r["sse"], ref["sse"]))
    if chk:
        same = sum(1 for c in chk if c[4])
        rel = max(abs(float(c[5]) - float(c[6])) / abs(float(c[6])) for c in chk if float(c[6]) != 0)
        L.append(f"{same}/{len(chk)} tuned/wilber cells have final SSE bit-identical to v8 "
                 f"(max relative diff {rel:.2e}).")
        for c in chk:
            if not c[4]:
                L.append(f"- diff {c[0]} {c[1]} k={c[2]} q={c[3] or '-'}: {c[5]} vs v8 {c[6]}")
    else:
        L.append("No matching v8 tuned/wilber cells (expected in --smoke on non-Mac platforms for normal/lognormal).")
    L += ["", f"Wall time: {wall/60:.1f} min. See runs.csv (one row per cell), run.log (raw driver output), ../env.txt."]
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
