#!/usr/bin/env python3
"""Orchestrate the LEX sig5 k-grid rerun on the reported (group-Refine) path. Invoked by run/kgrid_lex_fix.sh."""
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
OUT = PKG / "out" / "lex"
BUILD = PKG / "build"

# Fill cells only. Already-measured k in {64,256,512} are not re-run.
FILL_KS = (16, 32, 128)
SMALL_KS = (4, 8)
SMALL_K_DATASET_FULL = "normal_n1048576"
FULL_N = (1_048_576, 8_388_608)
SMOKE_N = (4_096,)
METHODS = ("sig5",)  # reported LEX path: 5-sig grouping, Refine on group atoms
FAMILIES = ("uniform", "normal", "lognormal")

AUDIO_GZ_NAME = "dharini-abhogi.txt.gz"
AUDIO_GZ_SHA = "3576f34420b749818b1897971692354fa6ded5af23a456ffba1965e1d8725c27"
AUDIO_TXT_SHA = "b19a45a308ba4758d8769b0d71e43ccf7503c0a0182e8d05e3ac5583586b2c3e"
AUDIO_N = 8_168_256

CSV_FIELDS = [
    "dataset", "method", "k", "M", "n", "D", "Dprime", "sig", "status",
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


def ks_for(name: str, smoke: bool) -> list[int]:
    if smoke:
        ks = [16, 32]
        if name.startswith("normal_"):
            return list(SMALL_KS) + ks
        return ks
    ks = list(FILL_KS)
    if name == SMALL_K_DATASET_FULL:
        return list(SMALL_KS) + ks
    return ks


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


def parse_exp_from_ns(ns: tuple[int, ...]) -> str:
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
        log(f"synth ready (smoke sizes {list(ns)}, hashes not checked against the shipped manifest)")
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


def write_summary(rows: list[dict], path: Path, smoke: bool) -> None:
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
    lines.append("# k-grid LEX fix: lex1dx sig5 (Refine on precision-group atoms)")
    lines.append("")
    if smoke:
        lines.append(
            "SMOKE: synthetic N = 4096, all three families, LEX sig5 at k ∈ {16,32} "
            "(plus k = 4, 8 on normal), one warm-up and one timed rep. "
            "Times are that single rep, not a median."
        )
    else:
        lines.append(
            "One warm-up and 3 timed reps; times are medians, in seconds. "
            "**Comparison is post-sort only**: every method is given the same array, "
            "already sorted with `std::sort` (that sort is timed separately and logged "
            "as `sort_s`; it is not added to `post_s`). "
            "lex1dx **sig5** always groups to 5 significant digits and Refines on "
            "the precision-group atoms; bounds are mapped to raw indices only for SSE "
            "scoring (the reported configuration). "
            "M = max(4096, 32k): M stays 4096 for k ≤ 128 and grows with k above that, "
            "so these cells evaluate the default budget policy rather than isolating k. "
            "If an exact SSE* is supplied it is the two-pass long-double SSE of unweighted "
            "Grønlund wilber-interp on the sorted points, and REE is (SSE − SSE*) / SSE* × 100. "
            "This pack does not rerun the exact solver. Peak RSS is the cell process (`getrusage`)."
        )
    lines.append("")
    lines.append("This pack reruns only LEX sig5. Exact-solver SSE* and REE are left blank unless an exact value was passed in.")
    lines.append("")
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

    def ds_key(d: str):
        if "_n" in d:
            try:
                return (0, int(d.split("_n")[-1]), d)
            except ValueError:
                return (1, 0, d)
        return (1, 0, d)

    datasets = sorted({k[0] for k in by}, key=ds_key)
    ks = sorted({k[2] for k in by})
    for ds in datasets:
        for k in ks:
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
    lines.append("See `out/lex/runs.csv` and `out/lex/run.log`. Re-run the same command to resume; finished (status=ok) cells are skipped.")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    smoke = "--smoke" in sys.argv
    python = os.environ.get("PYTHON_BIN") or sys.executable
    cell = find_bin("scale_cell")
    txt64 = find_bin("txt_to_f64")
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "runs.csv"
    log_path = OUT / "run.log"
    rows = load_csv(csv_path)

    ns = SMOKE_N if smoke else FULL_N
    reps = 1 if smoke else 3
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

    log(f"datasets={[d for d,_ in datasets]} reps={reps} smoke={int(smoke)}")
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

        for k in ks_for(name, smoke):
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
                sr = next((r for r in rows if r.get("dataset") == name and r.get("method") == "sort"), None)
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
    log(f"LEX sig5 fill done in {time.time()-t_all:.1f}s. Results: out/lex/summary.md")


if __name__ == "__main__":
    main()
