// One scaling cell, one process (so peak RSS is per cell).
//   scale_cell METHOD SORTED.f64 DATASET K REPS [EXACT_SSE] [--q=Q]
// METHOD: tuned | sig5 | sig5exact | wilber
//   sig5exact (v9 control) 5-sig grouping only (no k-fallback), weighted exact
//   wilber DP on all D' groups as point masses; no Collapse, no Refine.
//          as point masses; no Collapse, no Refine
//   --q=Q  Collapse quantile (default 0.75 = v8 behaviour)
// Built twice: scale_cell (patched bench/weighted_wilber) and, with -DLEX_PRISTINE,
// pristine_cell (unmodified vendor/kmeans1d; only method wilber, reported as
// wilber_pristine). Separate binaries, so the two kmeans_wilber never meet.
//   tuned  lex1dx tuned Collapse on full-precision atoms + weighted wilber DP + Refine
//   sig5   group to 5 significant digits, then the same pipeline (D' and sig logged).
//          If 5-sig D' < k and the raw data has D >= k, retry at 6, 7, ... digits
//          (up to 15), then raw distinct-value atoms, so k_eff = min(k, D).
//          After the DP, map bounds to raw indices and Refine on raw atoms.
//   wilber Grønlund wilber-interp on the raw sorted points
// Comparison is POST-SORT ONLY: input is sorted once (std::sort, timed separately as a
// reference column), every method gets the same sorted array. M = max(4096, 32k).
// One unrecorded warm-up + REPS timed runs; times are medians. Prints one CELL line.
#include "kmeans.h"
#ifndef LEX_PRISTINE
#include "lex/lex.hpp"
#include "detail.hpp"
#include "tuned_collapse.hpp"
#include "siggroup.hpp"
#define WILBER_NAME "wilber"
#else
#define WILBER_NAME "wilber_pristine"
#endif

#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using clk = std::chrono::steady_clock;
static double secs(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}
static double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
static double peak_rss_mb() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0); // bytes
#else
    return static_cast<double>(ru.ru_maxrss) / 1024.0;            // KiB
#endif
}

static std::vector<double> load_f64(const char* path) {
    FILE* fp = std::fopen(path, "rb");
    if (!fp) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(2); }
    std::fseek(fp, 0, SEEK_END);
    const long bytes = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);
    if (bytes <= 0 || bytes % 8 != 0) { std::fprintf(stderr, "bad f64 %s\n", path); std::exit(2); }
    const std::size_t n = static_cast<std::size_t>(bytes) / 8;
    std::vector<double> x(n);
    if (std::fread(x.data(), 8, n, fp) != n) { std::fprintf(stderr, "short read %s\n", path); std::exit(2); }
    std::fclose(fp);
    for (double v : x)
        if (!std::isfinite(v)) { std::fprintf(stderr, "non-finite value in %s\n", path); std::exit(2); }
    return x;
}

static double raw_sse(const std::vector<double>& sorted, const std::vector<int>& b) {
    long double sse = 0;
    for (std::size_t c = 0; c + 1 < b.size(); ++c) {
        const int lo = b[c], hi = b[c + 1];
        if (hi <= lo) continue;
        long double s = 0;
        for (int i = lo; i < hi; ++i) s += sorted[i];
        const long double m = s / (hi - lo);
        for (int i = lo; i < hi; ++i) {
            const long double d = sorted[i] - m;
            sse += d * d;
        }
    }
    return static_cast<double>(sse);
}

#ifndef LEX_PRISTINE
static std::vector<int> wilber_to_bounds(const std::vector<int>& groups,
                                         const std::vector<std::size_t>& path) {
    std::vector<std::size_t> p = path;
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    std::vector<int> b;
    const int G = static_cast<int>(groups.size()) - 1;
    for (std::size_t gi : p)
        if (gi <= static_cast<std::size_t>(G)) b.push_back(groups[static_cast<int>(gi)]);
    if (b.empty() || b.front() != 0) b.insert(b.begin(), 0);
    if (b.back() != groups[G]) b.push_back(groups[G]);
    return b;
}

struct Agg { std::vector<double> w, wx, wx2; };

static Agg agg_atoms(const lex::detail::Atoms& a, const std::vector<int>& groups) {
    Agg g;
    const int G = static_cast<int>(groups.size()) - 1;
    g.w.resize(G); g.wx.resize(G); g.wx2.resize(G);
    for (int i = 0; i < G; ++i) {
        const int lo = groups[i], hi = groups[i + 1];
        if (a.weighted) {
            // Point mass at the group mean: Wilber's interval cost drops a lone
            // item's internal SSE, so it must be excluded everywhere.
            const double W = a.wsum(lo, hi), mu = a.mean(lo, hi);
            g.w[i] = W; g.wx[i] = W * mu; g.wx2[i] = W * mu * mu;
        } else {
            const double c = hi - lo;
            const double mu = a.mean(lo, hi);
            g.w[i] = c; g.wx[i] = c * mu; g.wx2[i] = c * mu * mu;
        }
    }
    return g;
}

struct LexRun {
    double group_s = 0, atoms_s = 0, collapse_s = 0, dp_s = 0, refine_s = 0, post_s = 0;
    double sse = 0, coarse_sse = 0;
    int D = 0, Dprime = 0, sig = 0, atoms = 0, G = 0, k_eff = 0;
    int refine_sweeps = 0, refine_moves = 0;
    CollapseStats cstats;
    bool collapsed = false;
};

// `sorted` already sorted and -0 canonical.
// sig5 (library-matching): start at 5 significant digits; if D' < k and D >= k,
// raise the digit count (up to 15) or fall back to raw distinct-value atoms.
// After the DP, map group bounds to raw indices and Refine on raw atoms.
// sig5exact: always 5-sig groups, no fallback, no Collapse, no Refine.
static LexRun run_lex(const std::vector<double>& sorted, int k, std::size_t M, bool sig5, int D,
                      double q = 0.75, bool exact_only = false) {
    LexRun r;
    r.D = D;
    const auto t_all = clk::now();
    std::vector<int> elem_start; // group j covers sorted[elem_start[j] .. elem_start[j+1])
    lex::detail::Atoms* storage = nullptr;
    bool map_groups = false;

    if (sig5) {
        auto t0 = clk::now();
        int sig = 5;
        auto g = siggroup::group_fast(sorted, sig);
        bool use_raw = false;
        if (!exact_only) {
            while (static_cast<int>(g.w.size()) < k &&
                   static_cast<int>(g.w.size()) < D &&
                   sig < 15) {
                ++sig;
                g = siggroup::group_fast(sorted, sig);
            }
            if (static_cast<int>(g.w.size()) < k && static_cast<int>(g.w.size()) < D)
                use_raw = true;
        }
        r.group_s = secs(t0, clk::now());
        if (use_raw) {
            r.Dprime = D;
            r.sig = 0;
            t0 = clk::now();
            storage = new lex::detail::Atoms(sorted, 0.0);
            r.atoms_s = secs(t0, clk::now());
        } else {
            r.Dprime = static_cast<int>(g.w.size());
            r.sig = sig;
            elem_start = std::move(g.start);
            t0 = clk::now();
            // means are strictly increasing on sorted input, so sort_pairs is a no-op check
            lex::detail::sort_pairs(g.mean, g.w);
            storage = new lex::detail::Atoms(std::move(g.mean), std::move(g.w), 0.0);
            r.atoms_s = secs(t0, clk::now());
            map_groups = true;
        }
    } else {
        r.Dprime = D;
        r.sig = 0;
        auto t0 = clk::now();
        storage = new lex::detail::Atoms(sorted, 0.0);
        r.atoms_s = secs(t0, clk::now());
    }
    const lex::detail::Atoms& a = *storage;
    r.atoms = a.d;
    const int kk = std::min(k, a.d);
    r.k_eff = kk;
    const int m = static_cast<int>(std::max<std::size_t>(kk, std::min<std::size_t>(M, a.d)));

    auto t0 = clk::now();
    std::vector<int> groups;
    if (m >= a.d || exact_only) {
        groups = a.start;
    } else {
        groups = tuned_collapse(a, m, &r.cstats, q);
        r.collapsed = true;
    }
    r.collapse_s = secs(t0, clk::now());
    r.G = static_cast<int>(groups.size()) - 1;

    t0 = clk::now();
    std::vector<int> b;
    if (kk <= 1) {
        b = {0, a.n};
    } else if (kk >= r.G) {
        b = groups;
    } else {
        const Agg agg = agg_atoms(a, groups);
        kmeans_wilber fn(agg.w, agg.wx, agg.wx2);
        fn.set_search_strategy(search_strategy::INTERPOLATION);
        auto res = fn.compute_and_report(static_cast<std::size_t>(kk));
        if (!res) { std::fprintf(stderr, "weighted wilber DP failed\n"); std::exit(3); }
        b = wilber_to_bounds(groups, res->path);
    }
    r.dp_s = secs(t0, clk::now());

    if (map_groups) {
        std::vector<int> raw(b.size());
        for (std::size_t c = 0; c < b.size(); ++c) raw[c] = elem_start[b[c]];
        b = std::move(raw);
    }
    r.coarse_sse = raw_sse(sorted, b);

    if (!exact_only) {
        t0 = clk::now();
        if (map_groups) {
            lex::detail::Atoms raw_atoms(sorted, 0.0);
            const auto st = lex::detail::refine(raw_atoms, b, lex::Options{}.tolerance);
            r.refine_sweeps = st.sweeps;
            r.refine_moves = st.moves;
        } else {
            const auto st = lex::detail::refine(a, b, lex::Options{}.tolerance);
            r.refine_sweeps = st.sweeps;
            r.refine_moves = st.moves;
        }
        r.refine_s = secs(t0, clk::now());
    }
    r.post_s = secs(t_all, clk::now());
    r.sse = raw_sse(sorted, b);
    delete storage;
    return r;
}

#endif // LEX_PRISTINE

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: scale_cell sort IN.f64 OUT.f64 DATASET REPS\n"
            "       scale_cell {tuned|sig5|wilber} SORTED.f64 DATASET K REPS [EXACT_SSE]\n");
        return 2;
    }
    const std::string method = argv[1];

    if (method == "sort") {
        if (argc < 6) { std::fprintf(stderr, "usage: scale_cell sort IN.f64 OUT.f64 DATASET REPS\n"); return 2; }
        const char* path = argv[2];
        const char* outp = argv[3];
        const char* dataset = argv[4];
        const int reps = std::max(1, std::atoi(argv[5]));
        const std::vector<double> raw = load_f64(path);
        const std::size_t n = raw.size();
        std::vector<double> sort_t, sorted;
        for (int r = -1; r < reps; ++r) {
            std::vector<double> v = raw;
            const auto t0 = clk::now();
            for (double& x : v) x = x + 0.0;
            std::sort(v.begin(), v.end());
            const double t = secs(t0, clk::now());
            if (r >= 0) sort_t.push_back(t);
            if (r == -1) sorted = std::move(v);
        }
        int D = n ? 1 : 0;
        for (std::size_t i = 1; i < n; ++i) if (sorted[i] != sorted[i - 1]) ++D;
        FILE* out = std::fopen(outp, "wb");
        if (!out) { std::fprintf(stderr, "cannot write %s\n", outp); return 2; }
        if (std::fwrite(sorted.data(), 8, n, out) != n) return 2;
        std::fclose(out);
        std::printf("CELL dataset=%s method=sort k=0 M=0 n=%zu D=%d Dprime=%d status=ok "
                    "sort_s=%.6f post_s=0 peak_rss_mb=%.1f\n",
                    dataset, n, D, D, median(sort_t), peak_rss_mb());
        std::fflush(stdout);
        return 0;
    }

    if (argc < 6) {
        std::fprintf(stderr, "usage: scale_cell {tuned|sig5|wilber} SORTED.f64 DATASET K REPS [EXACT_SSE]\n");
        return 2;
    }
    const char* path = argv[2];
    const char* dataset = argv[3];
    const int k = std::atoi(argv[4]);
    const int reps = std::max(1, std::atoi(argv[5]));
    double exact_in = 0.0, q = 0.75;
    for (int i = 6; i < argc; ++i) {
        if (std::strncmp(argv[i], "--q=", 4) == 0) q = std::strtod(argv[i] + 4, nullptr);
        else exact_in = std::strtod(argv[i], nullptr);
    }
    if (!(q > 0.0 && q < 1.0)) { std::fprintf(stderr, "--q must be in (0,1)\n"); return 2; }
#ifdef LEX_PRISTINE
    if (method != "wilber") {
#else
    if (method != "tuned" && method != "sig5" && method != "sig5exact" && method != "wilber") {
#endif
        std::fprintf(stderr, "unknown method %s\n", method.c_str());
        return 2;
    }
    const std::size_t M = static_cast<std::size_t>(std::max(4096, 32 * k));

    std::vector<double> sorted = load_f64(path);
    const std::size_t n = sorted.size();
    // Caller hands a file already sorted by `scale_cell sort` (std::sort, -0 canonical).
    for (double& x : sorted) x = x + 0.0;
    if (!std::is_sorted(sorted.begin(), sorted.end())) {
        std::fprintf(stderr, "input is not sorted: %s\n", path);
        return 2;
    }
    int D = n ? 1 : 0;
    for (std::size_t i = 1; i < n; ++i) if (sorted[i] != sorted[i - 1]) ++D;
    const double sort_s = 0; // filled by the orchestrator from the sort cell

    if (method == "wilber") {
        if (static_cast<std::size_t>(k) >= n || k > D) {
            std::printf("CELL dataset=%s method=" WILBER_NAME " k=%d M=0 n=%zu D=%d status=skip_k_ge_D "
                        "sort_s=%.6f peak_rss_mb=%.1f\n", dataset, k, n, D, sort_s, peak_rss_mb());
            std::fflush(stdout);
            return 0;
        }
        std::vector<double> ts;
        double cost = 0, sse = 0;
        std::size_t path_n = 0;
        for (int r = -1; r < reps; ++r) {
            const auto t0 = clk::now();
            kmeans_wilber fn(sorted);
            fn.set_search_strategy(search_strategy::INTERPOLATION);
            auto res = fn.compute_and_report(static_cast<std::size_t>(k));
            const double t = secs(t0, clk::now());
            if (!res) {
                std::printf("CELL dataset=%s method=" WILBER_NAME " k=%d M=0 n=%zu D=%d status=fail_null "
                            "sort_s=%.6f peak_rss_mb=%.1f\n", dataset, k, n, D, sort_s, peak_rss_mb());
                return 4;
            }
            if (r >= 0) ts.push_back(t);
            if (r == reps - 1) {
                std::vector<int> b;
                for (auto it = res->path.rbegin(); it != res->path.rend(); ++it)
                    b.push_back(static_cast<int>(*it));
                cost = res->cost;
                sse = raw_sse(sorted, b);
                path_n = res->path.size();
            }
        }
        const int k_eff = static_cast<int>(path_n) - 1;
        std::printf("CELL dataset=%s method=" WILBER_NAME " k=%d M=0 n=%zu D=%d Dprime=%d status=ok "
                    "sort_s=%.6f post_s=%.6f dp_s=%.6f group_s=0 atoms_s=0 collapse_s=0 refine_s=0 "
                    "sse=%.17g lib_cost=%.17g exact_sse=%.17g ree_pct=0 coarse_sse=%.17g "
                    "k_eff=%d atoms=%zu G=%zu collapse_passes=0 collapse_merges=0 merges_per_pass=0 "
                    "refine_sweeps=0 refine_moves=0 reps=%d peak_rss_mb=%.1f\n",
                    dataset, k, n, D, D, sort_s, median(ts), median(ts),
                    sse, cost, sse, sse, k_eff, n, n, reps, peak_rss_mb());
        std::fflush(stdout);
        return 0;
    }

#ifdef LEX_PRISTINE
    (void)exact_in; (void)M;
    return 2;
#else
    const bool exact_only = (method == "sig5exact");
    const bool sig5 = (method == "sig5") || exact_only;
    (void)run_lex(sorted, k, M, sig5, D, q, exact_only); // warm-up
    std::vector<double> post, grp, atm, col, dp, ref;
    LexRun last;
    for (int r = 0; r < reps; ++r) {
        LexRun x = run_lex(sorted, k, M, sig5, D, q, exact_only);
        post.push_back(x.post_s); grp.push_back(x.group_s); atm.push_back(x.atoms_s);
        col.push_back(x.collapse_s); dp.push_back(x.dp_s); ref.push_back(x.refine_s);
        last = std::move(x);
    }
    const double ree = (exact_in > 0) ? (last.sse - exact_in) / exact_in * 100.0 : NAN;
    const double cree = (exact_in > 0) ? (last.coarse_sse - exact_in) / exact_in * 100.0 : NAN;
    std::printf("CELL dataset=%s method=%s k=%d M=%zu n=%zu D=%d Dprime=%d sig=%d status=ok "
                "sort_s=%.6f post_s=%.6f dp_s=%.6f group_s=%.6f atoms_s=%.6f collapse_s=%.6f refine_s=%.6f "
                "sse=%.17g lib_cost=0 exact_sse=%.17g ree_pct=%.6g coarse_sse=%.17g coarse_ree_pct=%.6g "
                "k_eff=%d atoms=%d G=%d collapse_passes=%d collapse_merges=%d merges_per_pass=%s "
                "refine_sweeps=%d refine_moves=%d reps=%d peak_rss_mb=%.1f q=%.4g\n",
                dataset, method.c_str(), k, M, n, D, last.Dprime, last.sig,
                sort_s, median(post), median(dp), median(grp), median(atm), median(col), median(ref),
                last.sse, exact_in, ree, last.coarse_sse, cree,
                last.k_eff, last.atoms, last.G, last.cstats.passes, last.cstats.merges_total(),
                last.cstats.merges_csv().c_str(),
                last.refine_sweeps, last.refine_moves, reps, peak_rss_mb(), q);
    std::fflush(stdout);
    return 0;
#endif
}
