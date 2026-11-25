"""v4 figures and tables for the 2026-11-25 LEX paper, full k grid.

Sources (read only; all paths are under this repository):
  Main N-grid at k in {64,256,512}: results/mac_out_v8/out/scaling/runs.csv
    methods: patched exact (wilber), full-precision LEX (tuned),
    precision-grouped LEX (sig5, group-level Refine).
  Fill-in k in {16,32,128} and small-k {4,8}:
    wilber from results/mac_out_kgrid_fill/out/lex/runs.csv;
    precision-grouped LEX (sig5, group-level Refine) from
    results/mac_out_kgrid_lex_fix/out/lex/runs.csv, with REE
    computed against the kgrid_fill wilber SSE*.
  Controls (stock exact, grouped exact, precision-grouped LEX) at k in
    {64,256,512} on the four large inputs:
    results/mac_out_v9/out/controls/runs.csv.
  Ckmeans.1d.dp 4.3.6, latest ok sorted row per cell, from
    results/ckmeans_matched, results/ckmeans_matched_k512,
    results/ckmeans_matched_k64_256, and results/mac_out_kgrid_fill/out/ckmeans.
  Variant timings (appendix): results/gronlund/mac_out_2026-11-24/out/results.csv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fractions import Fraction as F

# Package root (lex1dx-paper/). Paths stay inside this repository.
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "paper" / "figs"
TAB = ROOT / "paper" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

v8 = pd.read_csv(RES / "mac_out_v8/out/scaling/runs.csv")
v9 = pd.read_csv(RES / "mac_out_v9/out/controls/runs.csv")
kgrid_fill = pd.read_csv(RES / "mac_out_kgrid_fill/out/lex/runs.csv")
kgrid_fix = pd.read_csv(RES / "mac_out_kgrid_lex_fix/out/lex/runs.csv")

def _ree_vs_wilber(sig5_df, wilber_df):
    """Attach REE against existing wilber SSE* (wilber.sse)."""
    w = wilber_df[wilber_df.method == "wilber"].set_index(["dataset", "k"])["sse"]
    out = sig5_df.copy()
    star = [float(w.loc[(d, k)]) for d, k in zip(out.dataset, out.k)]
    out["exact_sse"] = star
    out["ree_pct"] = (out["sse"].astype(float) - out["exact_sse"]) / out["exact_sse"] * 100
    return out

sig5_fix = _ree_vs_wilber(kgrid_fix[kgrid_fix.method == "sig5"], kgrid_fill)
kgrid = pd.concat([
    kgrid_fill[kgrid_fill.method.isin(["wilber", "sort"])],
    sig5_fix,
], ignore_index=True)

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.4,
    "lines.markersize": 5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

COL = {"wilber": "black", "tuned": "#1f4e9b", "sig5": "#c44e16", "ckmeans": "#2a7a3a"}
NM = {
    "wilber": "exact (table-free)",
    "tuned": "full-precision LEX",
    "sig5": "precision-grouped LEX",
    "ckmeans": "Ckmeans.1d.dp",
}
MK = {"wilber": "s", "tuned": "o", "sig5": "D", "ckmeans": "v"}
FAMS = ["uniform", "normal", "lognormal"]
FAMLAB = {"uniform": "Uniform", "normal": "Normal", "lognormal": "Lognormal"}
LARGE = ["uniform_n8388608", "normal_n8388608", "lognormal_n8388608", "audio"]
LARGELAB = {
    "uniform_n8388608": r"Uniform, $N=2^{23}$",
    "normal_n8388608": r"Normal, $N=2^{23}$",
    "lognormal_n8388608": r"Lognormal, $N=2^{23}$",
    "audio": r"Audio, $N=8{,}168{,}256$",
}
M1 = ["uniform_n1048576", "normal_n1048576", "lognormal_n1048576"]
M1LAB = {
    "uniform_n1048576": r"Uniform, $N=2^{20}$",
    "normal_n1048576": r"Normal, $N=2^{20}$",
    "lognormal_n1048576": r"Lognormal, $N=2^{20}$",
}
KS = [16, 32, 64, 128, 256, 512]
DSLAB = {
    "uniform_n8388608": "Uniform",
    "normal_n8388608": "Normal",
    "lognormal_n8388608": "Lognormal",
    "audio": "Audio",
    "uniform_n1048576": "Uniform",
    "normal_n1048576": "Normal",
    "lognormal_n1048576": "Lognormal",
}

syn = v8[v8.dataset.str.contains("_n")].copy()
syn["fam"] = syn.dataset.str.split("_n").str[0]

# Fill-in cells for max-REE (sig5 only; no full-precision at these k)
kg_syn = kgrid[(kgrid.method == "sig5") & (kgrid.dataset.str.contains("_n"))
               & (kgrid.k.isin(KS))].copy()
kg_syn["fam"] = kg_syn.dataset.str.split("_n").str[0]


def xticks_N(ax):
    ticks = [2 ** i for i in (16, 18, 20, 23)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([rf"$2^{{{i}}}$" for i in (16, 18, 20, 23)])
    ax.set_xlim(2 ** 15.7, 2 ** 23.4)


def method_legend(fig, methods, extra=None, y=-0.02):
    hs = [Line2D([], [], color=COL[m], marker=MK[m], ls="-", label=NM[m]) for m in methods]
    if extra:
        hs.append(extra)
    fig.legend(handles=hs, loc="lower center", ncol=len(hs), frameon=False,
               bbox_to_anchor=(0.5, y), fontsize=8)


def fmt_s(x):
    x = float(x)
    if x >= 10:
        return f"{x:.2f}"
    if x >= 1:
        return f"{x:.3f}"
    return f"{x:.3f}"


def fmt_ree(x):
    x = float(x)
    if abs(x) < 5e-10:
        return "0"
    return f"{x:.5f}".rstrip("0").rstrip(".") if x < 0.001 else f"{x:.5f}".rstrip("0").rstrip(".")


# ---------------------------------------------------------------------------
# Ten-value pipeline figure (updated DP label)
# ---------------------------------------------------------------------------
def sse_parts(parts):
    t = F(0)
    for b in parts:
        m = F(sum(b), len(b))
        t += sum((F(v) - m) ** 2 for v in b)
    return t

xvals = [9, 12, 14, 17, 19, 22, 26, 27, 32, 36]
rows = [
    ("Collapse to five\ncoarse groups",
     [[9, 12, 14], [17, 19], [22], [26, 27], [32, 36]],
     None, "#555555"),
    ("DP, cuts between\nthose groups",
     [[9, 12, 14, 17, 19], [22, 26, 27], [32, 36]],
     sse_parts([[9, 12, 14, 17, 19], [22, 26, 27], [32, 36]]),
     "#2a7a3a"),
    ("Refine: 19 joins\nthe middle bin",
     [[9, 12, 14, 17], [19, 22, 26, 27], [32, 36]],
     sse_parts([[9, 12, 14, 17], [19, 22, 26, 27], [32, 36]]),
     "#1f4e9b"),
    ("Collapse to three bins,\nthen Refine",
     [[9, 12, 14], [17, 19, 22], [26, 27, 32, 36]],
     sse_parts([[9, 12, 14], [17, 19, 22], [26, 27, 32, 36]]),
     "#c44e16"),
]
assert rows[1][2] == F(424, 5) and rows[2][2] == 83 and rows[3][2] == F(1081, 12)

fig, ax = plt.subplots(figsize=(6.5, 4.85))
ysp = 1.18
for r, (name, parts, s, c) in enumerate(rows):
    y = -r * ysp
    for j, b in enumerate(parts):
        ax.add_patch(plt.Rectangle(
            (b[0] - 0.85, y - 0.26), b[-1] - b[0] + 1.7, 0.52,
            facecolor=c, alpha=0.12 + 0.10 * (j % 2), edgecolor=c, lw=1.3, zorder=1))
        ax.scatter(b, [y] * len(b), s=40, color=c, zorder=3, edgecolors="white", linewidths=0.3)
    ax.text(6.8, y + 0.30, name, fontsize=8.5, ha="left", va="bottom", color="k",
            zorder=5, linespacing=1.05)
    if s is None:
        lab = ""
    elif s == 83:
        lab = f"SSE {float(s):.1f}  (optimum)"
    elif s == F(1081, 12):
        lab = f"SSE {float(s):.2f}  (8.5% above)"
    else:
        lab = f"SSE {float(s):.1f}"
    if lab:
        ax.text(38.2, y, lab, fontsize=8.5, va="center", ha="left", color=c)

ax.plot([18, 18], [0.28, -3 * ysp - 0.26], color="k", ls=":", lw=1.2, zorder=4)
ax.text(18, 0.90, "optimal cut  17 | 19,  hidden inside a coarse group",
        fontsize=8, ha="center", va="bottom")
y_dp, y_ref = -1 * ysp, -2 * ysp
ax.annotate("", xy=(19, y_ref + 0.18), xytext=(19, y_dp - 0.18),
            arrowprops=dict(arrowstyle="->", color="k", lw=1.15), zorder=5)
ax.text(20.15, (y_dp + y_ref) / 2, r"net SSE change: $-1.8$", fontsize=8, va="center")
y_lab = -3 * ysp - 0.52
for v in xvals:
    dx = 0
    if v == 26:
        dx = -0.55
    elif v == 27:
        dx = 0.55
    ax.text(v + dx, y_lab, str(v), ha="center", va="top", fontsize=9)
ax.set_xlim(6.3, 49.0)
ax.set_ylim(y_lab - 0.45, 1.28)
ax.axis("off")
fig.tight_layout()
fig.savefig(FIG / "v4_pipeline.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Seven-number before/after (unchanged arithmetic)
# ---------------------------------------------------------------------------
seven = [3, 4, 6, 7, 9, 14, 22]
seven_rows = [
    ("Before",
     [[3, 4], [6, 7, 9], [14, 22]], F(223, 6), "#555555"),
    ("After the paired move",
     [[3, 4, 6], [7, 9, 14], [22]], F(92, 3), "#1f4e9b"),
]
assert seven_rows[0][2] == F(223, 6) and seven_rows[1][2] == F(92, 3)

fig, ax = plt.subplots(figsize=(3.25, 2.85))
ysp = 1.35
for r, (name, parts, s, c) in enumerate(seven_rows):
    y = -r * ysp
    for j, b in enumerate(parts):
        ax.add_patch(plt.Rectangle(
            (b[0] - 0.45, y - 0.28), b[-1] - b[0] + 0.9, 0.56,
            facecolor=c, alpha=0.12 + 0.10 * (j % 2), edgecolor=c, lw=1.3, zorder=1))
        ax.scatter(b, [y] * len(b), s=48, color=c, zorder=3, edgecolors="white", linewidths=0.3)
    ax.text(1.35, y, name, fontsize=7.6, ha="right", va="center")
    ax.text(24.6, y, rf"SSE ${s.numerator}/{s.denominator}$", fontsize=8, va="center", color=c)
ax.annotate("", xy=(6, -ysp + 0.18), xytext=(6, -0.18),
            arrowprops=dict(arrowstyle="->", color="k", lw=1.1), zorder=5)
ax.annotate("", xy=(14, -ysp + 0.18), xytext=(14, -0.18),
            arrowprops=dict(arrowstyle="->", color="k", lw=1.1), zorder=5)
ax.text(15.5, -ysp / 2, r"drop $13/2$", fontsize=7.6, va="center")
y_lab = -ysp - 0.55
for v in seven:
    dx = 0
    if v == 3:
        dx = -0.35
    elif v == 4:
        dx = 0.35
    elif v == 6:
        dx = -0.35
    elif v == 7:
        dx = 0.35
    ax.text(v + dx, y_lab, str(v), ha="center", va="top", fontsize=9)
ax.text(6, 0.48, "6", fontsize=8, ha="center", va="bottom")
ax.text(14, 0.48, "14", fontsize=8, ha="center", va="bottom")
ax.set_xlim(-2.3, 31.0)
ax.set_ylim(y_lab - 0.4, 0.85)
ax.axis("off")
fig.tight_layout()
fig.savefig(FIG / "v4_seven.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Max REE over available k (main accuracy figure)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.55), sharey=True)
for a, fam in zip(axes, FAMS):
    d = syn[syn.fam == fam]
    # full-precision: k in {64,256,512} only
    sub = d[d.method == "tuned"]
    g = sub.groupby("n", as_index=False)["ree_pct"].max().sort_values("n")
    a.plot(g.n, g.ree_pct, "-", color=COL["tuned"], marker=MK["tuned"], zorder=3)
    # precision-grouped: v8 plus fill-in k at N=2^20 and 2^23
    parts = [d[d.method == "sig5"][["n", "ree_pct"]],
             kg_syn[kg_syn.fam == fam][["n", "ree_pct"]]]
    sg = pd.concat(parts, ignore_index=True).groupby("n", as_index=False)["ree_pct"].max().sort_values("n")
    a.plot(sg.n, sg.ree_pct, "-", color=COL["sig5"], marker=MK["sig5"], zorder=3)
    a.axhline(0.01, color="0.35", lw=0.8, ls="--")
    a.set_xscale("log", base=2)
    xticks_N(a)
    a.set_ylim(0, 0.0105)
    a.set_title(FAMLAB[fam])
    a.set_xlabel("N")
    a.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
    a.set_axisbelow(True)
axes[0].set_ylabel("maximum REE over $k$ (%)")
axes[0].text(2 ** 16.15, 0.00955, "0.01% target", fontsize=7.5, color="0.3", va="top")
method_legend(fig, ["tuned", "sig5"], y=-0.02)
fig.tight_layout(rect=(0, 0.10, 1, 1))
fig.savefig(FIG / "v4_ree_max.pdf")
plt.close()

# ---------------------------------------------------------------------------
# REE detail (appendix): full curves for 64,256,512; markers for 16,32,128
# ---------------------------------------------------------------------------
ls_k = {16: (0, (1, 2)), 32: (0, (2, 2)), 64: ":", 128: (0, (4, 2)), 256: "--", 512: "-"}
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.55), sharey=True)
for a, fam in zip(axes, FAMS):
    d = syn[syn.fam == fam]
    for m in ("tuned", "sig5"):
        for k in (64, 256, 512):
            s = d[(d.method == m) & (d.k == k)].sort_values("n")
            a.plot(s.n, s.ree_pct, ls_k[k], color=COL[m], marker=MK[m], markersize=4)
    # fill-in k as markers only (sig5)
    kgf = kg_syn[kg_syn.fam == fam]
    for k in (16, 32, 128):
        s = kgf[kgf.k == k]
        if len(s):
            a.plot(s.n, s.ree_pct, linestyle="None", color=COL["sig5"],
                   marker="x", markersize=5, markeredgewidth=1.0, zorder=4)
    a.axhline(0.01, color="0.35", lw=0.8, ls="--")
    a.set_xscale("log", base=2)
    a.set_yscale("log")
    xticks_N(a)
    a.set_ylim(2e-5, 2e-2)
    a.set_title(FAMLAB[fam])
    a.set_xlabel("N")
axes[0].set_ylabel("REE (%)")
h1 = [Line2D([], [], color=COL[m], marker=MK[m], label=NM[m]) for m in ("tuned", "sig5")]
h2 = [Line2D([], [], color="0.3", ls=ls_k[k], label=f"$k = {k}$") for k in (64, 256, 512)]
h2.append(Line2D([], [], color=COL["sig5"], marker="x", ls="None",
                 label=r"$k=16,32,128$ (precision-grouped)"))
fig.legend(handles=h1, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.08), fontsize=8)
fig.legend(handles=h2, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02), fontsize=7.5)
fig.tight_layout(rect=(0, 0.16, 1, 1))
fig.savefig(FIG / "v4_ree_detail.pdf")
plt.close()

# ---------------------------------------------------------------------------
# 2x3 time and RSS at k = 512 (hide inner y tick labels)
# ---------------------------------------------------------------------------
K = 512
fig, axes = plt.subplots(2, 3, figsize=(6.5, 7.0), sharex="col")
for j, fam in enumerate(FAMS):
    d = syn[(syn.fam == fam) & (syn.k.isin([0, K]))]
    ax_t, ax_r = axes[0, j], axes[1, j]
    for m in ("wilber", "tuned", "sig5"):
        s = d[d.method == m].sort_values("n")
        ax_t.plot(s.n, s.post_s, "-", color=COL[m], marker=MK[m])
        ax_r.plot(s.n, s.peak_rss_mb, "-", color=COL[m], marker=MK[m])
    srt = d[d.method == "sort"].sort_values("n")
    ax_t.plot(srt.n, srt.sort_s, "-", color="0.55", lw=2.0, alpha=0.85, zorder=2)
    ax_t.set_xscale("log", base=2)
    ax_t.set_yscale("log")
    ax_r.set_xscale("log", base=2)
    ax_r.set_yscale("log")
    xticks_N(ax_t)
    xticks_N(ax_r)
    ax_t.set_title(FAMLAB[fam])
    ax_r.set_xlabel("N")
    ax_t.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
    ax_r.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
    ax_t.set_axisbelow(True)
    ax_r.set_axisbelow(True)
    if j > 0:
        ax_t.tick_params(labelleft=False)
        ax_r.tick_params(labelleft=False)
axes[0, 0].set_ylabel("time after sorting (s)")
axes[1, 0].set_ylabel("peak RSS (MiB)")
tmin = min(ax.get_ylim()[0] for ax in axes[0])
tmax = max(ax.get_ylim()[1] for ax in axes[0])
rmin = min(ax.get_ylim()[0] for ax in axes[1])
rmax = max(ax.get_ylim()[1] for ax in axes[1])
for ax in axes[0]:
    ax.set_ylim(tmin, tmax)
for ax in axes[1]:
    ax.set_ylim(rmin, rmax)
fig.suptitle(r"$k = 512$", fontsize=11, y=0.995)
method_legend(fig, ["wilber", "tuned", "sig5"],
              extra=Line2D([], [], color="0.55", lw=2.0, label="std::sort"),
              y=-0.01)
fig.tight_layout(rect=(0, 0.06, 1, 0.98))
fig.savefig(FIG / "v4_time_rss_k512.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Speedup vs N (appendix)
# ---------------------------------------------------------------------------
w = v8[v8.method == "wilber"].set_index(["dataset", "k"]).post_s
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.4), sharey=True)
for a, fam in zip(axes, FAMS):
    d = syn[syn.fam == fam]
    for m in ("tuned", "sig5"):
        for k in (64, 256, 512):
            s = d[(d.method == m) & (d.k == k)].sort_values("n")
            a.plot(s.n, [w[(r.dataset, r.k)] / r.post_s for r in s.itertuples()],
                   ls_k[k], color=COL[m], marker=MK[m], markersize=4)
    # fill-in k speedup markers
    kgf = kgrid[(kgrid.dataset.str.startswith(fam + "_n")) & (kgrid.method.isin(["sig5", "wilber"]))]
    ww = kgf[kgf.method == "wilber"].set_index(["dataset", "k"]).post_s
    for k in (16, 32, 128):
        s = kgf[(kgf.method == "sig5") & (kgf.k == k)]
        if len(s):
            ys = [ww[(r.dataset, r.k)] / r.post_s for r in s.itertuples()]
            a.plot(s.n, ys, linestyle="None", color=COL["sig5"],
                   marker="x", markersize=5, markeredgewidth=1.0, zorder=4)
    a.set_xscale("log", base=2)
    a.set_yscale("log")
    xticks_N(a)
    a.set_title(FAMLAB[fam])
    a.set_xlabel("N")
    a.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
    a.set_axisbelow(True)
axes[0].set_ylabel("speedup after sorting")
h1 = [Line2D([], [], color=COL[m], marker=MK[m], label=NM[m]) for m in ("tuned", "sig5")]
h2 = [Line2D([], [], color="0.3", ls=ls_k[k], label=f"$k = {k}$") for k in (64, 256, 512)]
h2.append(Line2D([], [], color=COL["sig5"], marker="x", ls="None",
                 label=r"$k=16,32,128$ (precision-grouped)"))
fig.legend(handles=h1, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.08), fontsize=8)
fig.legend(handles=h2, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02), fontsize=7.5)
fig.tight_layout(rect=(0, 0.16, 1, 1))
fig.savefig(FIG / "v4_speedup_vs_N.pdf")
plt.close()

# ---------------------------------------------------------------------------
# D'/D with keyed real-file markers
# ---------------------------------------------------------------------------
s5 = v8[v8.method == "sig5"].drop_duplicates("dataset")
fig, ax = plt.subplots(figsize=(6.5, 3.35))
for fam, mk in zip(FAMS, "os^"):
    s = s5[s5.dataset.str.startswith(fam + "_n")].sort_values("n")
    ax.plot(s.n, 100 * s.Dprime / s.D, "-", marker=mk, label=FAMLAB[fam], color="0.15",
            markerfacecolor="white", markeredgecolor="0.15", markersize=6)
real = [
    ("audio", "audio", "P", "#1f4e9b"),
    ("h17_c_attn", "GPT-2 h17", "X", "#c44e16"),
    ("h35_c_proj", "GPT-2 h35", "+", "#2a7a3a"),
]
for ds, lab, mk, c in real:
    r = s5[s5.dataset == ds].iloc[0]
    pct = 100 * r.Dprime / r.D
    ax.plot(r.n, pct, marker=mk, color=c, markersize=8, markeredgewidth=1.4,
            linestyle="None", label=f"{lab} ({pct:.1f}%)", zorder=5)
ax.set_xscale("log", base=2)
ax.set_yscale("log")
xticks_N(ax)
yt = [100, 50, 20, 10, 5, 3]
ax.set_yticks(yt)
ax.set_yticklabels([str(t) for t in yt])
ax.set_ylim(2.4, 120)
ax.set_xlabel("N")
ax.set_ylabel(r"$D_5/D$ (%)")
ax.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig(FIG / "v4_dprime.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Time vs k at N ≈ 8M (companion to Table 4)
# ---------------------------------------------------------------------------
def load_ckmeans():
    srcs = [
        RES / "ckmeans_matched/runs.csv",
        RES / "ckmeans_matched_k512/out/ckmeans_matched/runs.csv",
        RES / "ckmeans_matched_k64_256/out/ckmeans_matched/runs.csv",
        RES / "mac_out_kgrid_fill/out/ckmeans/runs.csv",
    ]
    ck = pd.concat([pd.read_csv(p) for p in srcs], ignore_index=True)
    ok = ck[(ck.method == "Ckmeans.1d.dp") & (ck.status == "ok") & (ck.input == "sorted")]
    return ok.groupby(["dataset", "k"], as_index=False).tail(1)

ck_last = load_ckmeans()


def lex_wilber_row(ds, k, method):
    if k in (64, 256, 512):
        src = v8
    else:
        src = kgrid
    hit = src[(src.dataset == ds) & (src.k == k) & (src.method == method)]
    if len(hit) != 1:
        raise SystemExit(f"expected 1 row for {ds} {method} k={k}, got {len(hit)}")
    return hit.iloc[0]


fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.6), sharex=True)
for ax, ds in zip(axes.ravel(), LARGE):
    csub = ck_last[ck_last.dataset == ds].set_index("k")
    ck_t = [csub.loc[k, "post_s"] for k in KS]
    lex_t = [lex_wilber_row(ds, k, "sig5").post_s for k in KS]
    w_t = [lex_wilber_row(ds, k, "wilber").post_s for k in KS]
    ax.plot(KS, ck_t, "-", color=COL["ckmeans"], marker=MK["ckmeans"], label=NM["ckmeans"])
    ax.plot(KS, w_t, "-", color=COL["wilber"], marker=MK["wilber"], label=NM["wilber"])
    ax.plot(KS, lex_t, "-", color=COL["sig5"], marker=MK["sig5"], label=NM["sig5"])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(KS)
    ax.set_xticklabels([str(k) for k in KS])
    ax.set_title(LARGELAB[ds], fontsize=9)
    ax.yaxis.grid(True, ls=":", lw=0.5, color="0.85")
    ax.set_axisbelow(True)
axes[0, 0].set_ylabel("time after sorting (s)")
axes[1, 0].set_ylabel("time after sorting (s)")
axes[1, 0].set_xlabel("k")
axes[1, 1].set_xlabel("k")
handles = [Line2D([], [], color=COL[m], marker=MK[m], ls="-", label=NM[m])
           for m in ("ckmeans", "wilber", "sig5")]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02), fontsize=8)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(FIG / "v4_time_vs_k.pdf")
plt.close()

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def v9_sig5(ds, k):
    sB = v9[(v9.dataset == ds) & (v9.k == k) & (v9.method == "sig5") & (v9.experiment == "B")]
    return sB.iloc[0]


def v9_row(ds, k, method):
    return v9[(v9.dataset == ds) & (v9.k == k) & (v9.method == method)].iloc[0]


# Table 3: methods in rows, grouped by dataset, k=512, four large inputs.
# Exact = unmodified (stock) solver from controls.
lines = [
    r"\setlength{\tabcolsep}{4pt}",
    r"\begin{tabular}{@{}llrrr@{}}",
    r"\toprule",
    r"Data & Method & Time (s) & Peak RSS (MiB) & REE (\%) \\",
    r"\midrule",
]
for ds in LARGE:
    stock = v9_row(ds, 512, "wilber_pristine")
    se = v9_row(ds, 512, "sig5exact")
    tu = v8[(v8.dataset == ds) & (v8.k == 512) & (v8.method == "tuned")].iloc[0]
    sg = v9_sig5(ds, 512)
    block = [
        ("Exact reference", stock.post_s, stock.peak_rss_mb, None, stock.post_s / stock.post_s),
        ("Exact solve on precision groups", se.post_s, se.peak_rss_mb, se.ree_pct, stock.post_s / se.post_s),
        ("Full-precision LEX", tu.post_s, tu.peak_rss_mb, tu.ree_pct, stock.post_s / tu.post_s),
        ("Precision-grouped LEX", sg.post_s, sg.peak_rss_mb, sg.ree_pct, stock.post_s / sg.post_s),
    ]
    first = True
    for name, t, rss, ree, sp in block:
        lab = DSLAB[ds] if first else ""
        first = False
        ree_s = "Reference" if ree is None else f"{ree:.5f}"
        sp_s = "" if name == "Exact reference" else f" ({sp:.1f}$\\times$)"
        # Exact-reference times in this compact table use two decimals,
        # matching the r3 manuscript (9.18 s, 8.86 s, 11.96 s, 11.60 s).
        # The appendix grid table keeps three decimals when t < 10.
        if name == "Exact reference":
            tfmt = f"{t:.2f}"
        else:
            tfmt = f"{t:.3f}" if t < 10 else f"{t:.2f}"
            tfmt = f"{tfmt}{sp_s}"
        lines.append(f"{lab} & {name} & {tfmt} & {rss:.0f} & {ree_s} \\\\")
    if ds != LARGE[-1]:
        lines.append(r"\addlinespace")
# GPT-2 rows: only the weighted table-free solver was run on these (v8), so it is their baseline.
lines.append(r"\midrule")
for ds in ("h17_c_attn", "h35_c_proj"):
    w = v8[(v8.dataset == ds) & (v8.k == 512) & (v8.method == "wilber")].iloc[0]
    tu = v8[(v8.dataset == ds) & (v8.k == 512) & (v8.method == "tuned")].iloc[0]
    sg = v8[(v8.dataset == ds) & (v8.k == 512) & (v8.method == "sig5")].iloc[0]
    lab = {"h17_c_attn": r"GPT-2 h17 $c_{\mathrm{attn}}$", "h35_c_proj": r"GPT-2 h35 $c_{\mathrm{proj}}$"}[ds]
    lines.append(f"{lab} & Exact reference$^\\dagger$ & {w.post_s:.2f} & {w.peak_rss_mb:.0f} & Reference \\\\")
    lines.append(f" & Full-precision LEX & {tu.post_s:.3f} ({w.post_s/tu.post_s:.1f}$\\times$) & {tu.peak_rss_mb:.0f} & {tu.ree_pct:.5f} \\\\")
    lines.append(f" & Precision-grouped LEX & {sg.post_s:.3f} ({w.post_s/sg.post_s:.1f}$\\times$) & {sg.peak_rss_mb:.0f} & {sg.ree_pct:.5f} \\\\")
    if ds != "h35_c_proj":
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_headline_k512.tex").write_text("\n".join(lines) + "\n")

# Appendix: full 3-k headline for four large inputs (stock exact, available methods)
# Keep k=64,256,512 here (tuned and sig5exact only exist at those k).
lines = [
    r"\setlength{\tabcolsep}{3.5pt}",
    r"\resizebox{\textwidth}{!}{%",
    r"\begin{tabular}{@{}ll rrr rrr rrr@{}}",
    r"\toprule",
    r" &  & \multicolumn{3}{c}{$k=64$} & \multicolumn{3}{c}{$k=256$} & \multicolumn{3}{c}{$k=512$} \\",
    r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
    r"Data & Method & s & RSS & REE & s & RSS & REE & s & RSS & REE \\",
    r"\midrule",
]
methods_app = [
    ("Exact reference", "wilber_pristine", "v9"),
    ("Exact solve on precision groups", "sig5exact", "v9"),
    ("Full-precision LEX", "tuned", "v8"),
    ("Precision-grouped LEX", "sig5", "v9"),
]
for ds in LARGE:
    first = True
    for name, meth, src in methods_app:
        cells = []
        for k in (64, 256, 512):
            if src == "v8":
                r = v8[(v8.dataset == ds) & (v8.k == k) & (v8.method == meth)].iloc[0]
            elif meth == "sig5":
                r = v9_sig5(ds, k)
            else:
                r = v9_row(ds, k, meth)
            stock = v9_row(ds, k, "wilber_pristine")
            t = r.post_s
            tfmt = f"{t:.3f}" if t < 10 else f"{t:.2f}"
            if meth != "wilber_pristine":
                tfmt += f" ({stock.post_s/t:.1f}$\\times$)"
            ree = "—" if meth == "wilber_pristine" else f"{r.ree_pct:.5f}"
            cells.append(f"{tfmt} & {r.peak_rss_mb:.0f} & {ree}")
        lab = DSLAB[ds] if first else ""
        first = False
        lines.append(f"{lab} & {name} & " + " & ".join(cells) + r" \\")
    if ds != LARGE[-1]:
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}%", r"}"]
(TAB / "v4_headline_full.tex").write_text("\n".join(lines) + "\n")

# Table 4: Ckmeans / LEX / exact on four large inputs, all 6 k. Memory in GiB.
def gib(mb):
    return float(mb) / 1024.0

lines = [
    r"\setlength{\tabcolsep}{5.5pt}",
    r"\begin{tabular}{@{}ll rr rrr rr@{}}",
    r"\toprule",
    r" &  & \multicolumn{2}{c}{Ckmeans.1d.dp} & \multicolumn{3}{c}{Precision-grouped LEX} & \multicolumn{2}{c}{Table-free exact} \\",
    r"\cmidrule(lr){3-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
    r"Data & $k$ & Time (s) & RSS (GiB) & Time (s) & RSS (GiB) & REE (\%) & Time (s) & RSS (GiB) \\",
    r"\midrule",
]
for ds in LARGE:
    first = True
    for k in KS:
        c = ck_last[(ck_last.dataset == ds) & (ck_last.k == k)].iloc[0]
        s = lex_wilber_row(ds, k, "sig5")
        w = lex_wilber_row(ds, k, "wilber")
        lab = DSLAB[ds] if first else ""
        first = False
        lines.append(
            f"{lab} & {k} & "
            f"{c.post_s:.2f} & {gib(c.peak_rss_mb):.2f} & "
            f"{s.post_s:.3f} & {gib(s.peak_rss_mb):.3f} & {s.ree_pct:.5f} & "
            f"{w.post_s:.2f} & {gib(w.peak_rss_mb):.3f} \\\\"
        )
    if ds != LARGE[-1]:
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_ckmeans_8M.tex").write_text("\n".join(lines) + "\n")

# Appendix Table 4 companion: N=2^20
lines = [
    r"\setlength{\tabcolsep}{5.5pt}",
    r"\begin{tabular}{@{}ll rr rrr rr@{}}",
    r"\toprule",
    r" &  & \multicolumn{2}{c}{Ckmeans.1d.dp} & \multicolumn{3}{c}{Precision-grouped LEX} & \multicolumn{2}{c}{Table-free exact} \\",
    r"\cmidrule(lr){3-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
    r"Data & $k$ & Time (s) & RSS (GiB) & Time (s) & RSS (GiB) & REE (\%) & Time (s) & RSS (GiB) \\",
    r"\midrule",
]
for ds in M1:
    first = True
    for k in KS:
        c = ck_last[(ck_last.dataset == ds) & (ck_last.k == k)].iloc[0]
        s = lex_wilber_row(ds, k, "sig5")
        w = lex_wilber_row(ds, k, "wilber")
        lab = DSLAB[ds] if first else ""
        first = False
        lines.append(
            f"{lab} & {k} & "
            f"{c.post_s:.3f} & {gib(c.peak_rss_mb):.2f} & "
            f"{s.post_s:.3f} & {gib(s.peak_rss_mb):.3f} & {s.ree_pct:.5f} & "
            f"{w.post_s:.3f} & {gib(w.peak_rss_mb):.3f} \\\\"
        )
    if ds != M1[-1]:
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_ckmeans_1M.tex").write_text("\n".join(lines) + "\n")

# REE summary (update counts for sig5: 81+21=102, plus note)
# Reuse family breakdown from v8 (81) and mention fill-in separately in the caption.
def fam_of(ds):
    if ds.startswith("uniform"):
        return "Uniform"
    if ds.startswith("normal"):
        return "Normal"
    if ds.startswith("lognormal"):
        return "Lognormal"
    return "Real"


def ree_block(df):
    out = {}
    for fam, g in df.groupby(df.dataset.map(fam_of)):
        out[fam] = (len(g), g.ree_pct.min(), g.ree_pct.median(), g.ree_pct.max())
    return out


sig5_all = pd.concat([
    v8[v8.method == "sig5"],
    kgrid[(kgrid.method == "sig5") & (kgrid.k.isin(KS))],
], ignore_index=True)
tuned_all = v8[v8.method == "tuned"]
# real files: audio in both, gpt2 only v8
lines = [
    r"\begin{tabular}{@{}l rrrr rrrr@{}}",
    r"\toprule",
    r" & \multicolumn{4}{c}{Full-precision LEX} & \multicolumn{4}{c}{Precision-grouped LEX} \\",
    r"\cmidrule(lr){2-5}\cmidrule(lr){6-9}",
    r"Family & $n$ & Min & Median & Max & $n$ & Min & Median & Max \\",
    r"\midrule",
]
# Order
order = ["Uniform", "Normal", "Lognormal", "Real"]
# Real = audio + gpt2 from v8, plus kgrid audio extra k
for fam in order:
    if fam != "Real":
        t = tuned_all[tuned_all.dataset.map(fam_of) == fam]
        s = sig5_all[sig5_all.dataset.map(fam_of) == fam]
    else:
        t = tuned_all[~tuned_all.dataset.str.contains("_n")]
        s = sig5_all[~sig5_all.dataset.str.contains("_n")]
    lines.append(
        f"{fam} & {len(t)} & {t.ree_pct.min():.5f} & {t.ree_pct.median():.5f} & {t.ree_pct.max():.5f} & "
        f"{len(s)} & {s.ree_pct.min():.5f} & {s.ree_pct.median():.5f} & {s.ree_pct.max():.5f} \\\\"
    )
t, s = tuned_all, sig5_all
lines.append(r"\midrule")
lines.append(
    f"All & {len(t)} & {t.ree_pct.min():.5f} & {t.ree_pct.median():.5f} & {t.ree_pct.max():.5f} & "
    f"{len(s)} & {s.ree_pct.min():.5f} & {s.ree_pct.median():.5f} & {s.ree_pct.max():.5f} \\\\"
)
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_ree_summary.tex").write_text("\n".join(lines) + "\n")

# Quantile and GPT-2: same sources as before, regenerated without run tags in the tex captions.
# Quantile from v9 experiment C
qrows = v9[v9.experiment == "C"].copy()
# Table already exists conceptually; rebuild
# Look at unique q, methods, datasets
lines = [
    r"\begin{tabular}{@{}ll rrr rrr@{}}",
    r"\toprule",
    r" &  & \multicolumn{3}{c}{Full-precision LEX} & \multicolumn{3}{c}{Precision-grouped LEX} \\",
    r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
    r"Data & $q$ & Time (s) & RSS (MiB) & REE (\%) & Time (s) & RSS (MiB) & REE (\%) \\",
    r"\midrule",
]
for ds in ("normal_n8388608", "audio"):
    first = True
    for q in (0.5, 0.75, 0.9):
        tu = qrows[(qrows.dataset == ds) & (qrows.method == "tuned") & (qrows.q == q)]
        sg = qrows[(qrows.dataset == ds) & (qrows.method == "sig5") & (qrows.q == q)]
        if not len(tu) or not len(sg):
            continue
        tu, sg = tu.iloc[0], sg.iloc[0]
        lab = "Normal" if ds.startswith("normal") else "Audio"
        lab = lab if first else ""
        first = False
        lines.append(
            f"{lab} & {q} & {tu.post_s:.3f} & {tu.peak_rss_mb:.0f} & {tu.ree_pct:.5f} & "
            f"{sg.post_s:.3f} & {sg.peak_rss_mb:.0f} & {sg.ree_pct:.5f} \\\\"
        )
    if ds != "audio":
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_quantile.tex").write_text("\n".join(lines) + "\n")

# GPT-2, v8, patched exact as baseline
lines = [
    r"\setlength{\tabcolsep}{3.5pt}",
    r"\begin{tabular}{@{}ll rrr rrr rrr@{}}",
    r"\toprule",
    r" &  & \multicolumn{3}{c}{Exact reference} & \multicolumn{3}{c}{Full-precision LEX} & \multicolumn{3}{c}{Precision-grouped LEX} \\",
    r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
    r"Matrix & $k$ & Time (s) & RSS &  & Time (s) & RSS & REE (\%) & Time (s) & RSS & REE (\%) \\",
    r"\midrule",
]
mlab = {"h17_c_attn": r"h17 $c_{\mathrm{attn}}$", "h35_c_proj": r"h35 $c_{\mathrm{proj}}$"}
for ds in ("h17_c_attn", "h35_c_proj"):
    first = True
    for k in (64, 256, 512):
        w = v8[(v8.dataset == ds) & (v8.k == k) & (v8.method == "wilber")].iloc[0]
        tu = v8[(v8.dataset == ds) & (v8.k == k) & (v8.method == "tuned")].iloc[0]
        sg = v8[(v8.dataset == ds) & (v8.k == k) & (v8.method == "sig5")].iloc[0]
        lab = mlab[ds] if first else ""
        first = False
        lines.append(
            f"{lab} & {k} & {w.post_s:.3f} & {w.peak_rss_mb:.0f} &  & "
            f"{tu.post_s:.3f} ({w.post_s/tu.post_s:.1f}$\\times$) & {tu.peak_rss_mb:.0f} & {tu.ree_pct:.5f} & "
            f"{sg.post_s:.3f} ({w.post_s/sg.post_s:.1f}$\\times$) & {sg.peak_rss_mb:.0f} & {sg.ree_pct:.5f} \\\\"
        )
    if ds != "h35_c_proj":
        lines.append(r"\addlinespace")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_gpt2.tex").write_text("\n".join(lines) + "\n")

# Appendix: kmeans1d variant driver times, audio k=512 and a compact grid.
g = pd.read_csv(RES / "gronlund/mac_out_2026-11-24/out/results.csv")
# Use in-process driver_s when present (Grønlund methods); LEX rows have only elapsed_s.
# For the appendix we only need the four Grønlund methods on audio at k=64,256,512
# plus a note on ties failures. Also include wilber-interp driver_s to match the paper's earlier 12.3 s.
var_names = {
    "wilber-interp": r"\code{wilber-interp}",
    "hirc": r"\code{hirc}",
    "wilber-binary": r"\code{wilber-binary}",
    "monotone:compute": r"\code{monotone}",
}
lines = [
    r"\begin{tabular}{@{}l rrr@{}}",
    r"\toprule",
    r"Method & $k=64$ & $k=256$ & $k=512$ \\",
    r"\midrule",
]
for meth, lab in var_names.items():
    cells = []
    for k in (64, 256, 512):
        r = g[(g.dataset == "audio") & (g.k == k) & (g.method == meth)].iloc[0]
        cells.append(f"{r.driver_s:.1f}")
    lines.append(f"{lab} & " + " & ".join(cells) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]
(TAB / "v4_variants_audio.tex").write_text("\n".join(lines) + "\n")

print("wrote figures and tables")
for p in sorted(FIG.glob("v4_*")):
    print(" ", p.name, p.stat().st_size)
for p in sorted(TAB.glob("v4_*")):
    print(" ", p.name, p.stat().st_size)
