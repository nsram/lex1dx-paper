"""Reproduce the time-vs-N figure from a benchmark results CSV.

    python3 plot_timing.py results.csv [-o timing.png] [--post-sort] [--section timing]

The CSV is the `results.csv` written by the benchmark pack (v7 format: one row per section, dist,
N, seed, k, method, rule, with time_s, post_sort_s and presorted_s in seconds). The benchmark
harness itself is not part of this repository. Lines are medians over the cells (distributions x
seeds) at each N; exact solvers are dashed. With --post-sort, LEX tiers and ckmeans_exact show their
time after the sort and the baselines their time on presorted input. Needs pandas and matplotlib.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

# (label, method and rule as written in the CSV, exact solver?). The CSV predates the v6 names:
# LEX_default_coarse = lex_hybrid, LEX_default = lex_plain, fit_exact = ckmeans_exact.
SERIES = [
    ("lex_hybrid (coarse 16384)", "LEX_default_coarse", "16384", False),
    ("lex_plain", "LEX_default", "off", False),
    ("lex_plain, deck (1e-4 grid)", "LEX_deck", "off", False),
    ("ckmeans_exact (DP)", "fit_exact", "none", True),
    ("Ckmeans.1d.dp (R)", "Ckmeans.1d.dp", "none", True),
    ("ckmeans (Python/Rust; k=255 for 256)", "ckmeans_py*", "none", True),
]
BASELINES = {"Ckmeans.1d.dp", "ckmeans_py", "ckmeans_py_k255"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="timing.png")
    ap.add_argument("--post-sort", action="store_true")
    ap.add_argument("--section", default="timing")
    args = ap.parse_args()

    d = pd.read_csv(args.csv)
    d = d[(d.section == args.section) & (d.start != "error")]
    # keep the k values the LEX tiers ran at (drops the extra k = 255 rows run for the Python ckmeans)
    d = d[d.k.isin(d[d.method.str.startswith("LEX")].k.unique())]
    d["rule"] = d["rule"].astype(str)
    ks = sorted(int(k) for k in d.k.unique())
    fig, axes = plt.subplots(1, len(ks), figsize=(4.2 * len(ks), 3.8), squeeze=False)
    rows = []
    for ax, k in zip(axes[0], ks):
        dk = d[d.k == k]
        for label, method, rule, exact in SERIES:
            m = dk[dk.method.str.startswith(method[:-1])] if method.endswith("*") else dk[dk.method == method]
            m = m[m.rule == rule]
            if m.empty:
                continue
            base = m.method.iloc[0].startswith(("Ckmeans", "ckmeans"))
            col = ("presorted_s" if base else "post_sort_s") if args.post_sort else "time_s"
            g = m.groupby("N")[col].median()
            ax.plot(g.index, g.values, "--o" if exact else "-o", ms=3, label=label)
            rows += [(k, label, int(n), t) for n, t in g.items()]
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"k = {k}"); ax.set_xlabel("N"); ax.grid(True, which="both", alpha=0.3)
    axes[0][0].set_ylabel("median seconds (" + ("after the sort / presorted" if args.post_sort else "total, unsorted input") + ")")
    axes[0][-1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    t = pd.DataFrame(rows, columns=["k", "method", "N", "median_s"]).pivot_table(index=["k", "N"], columns="method", values="median_s")
    print(t.to_string(float_format=lambda x: f"{x:.4g}"))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
