// Score a 1-d partition with the same two-pass long-double SSE as v10 scale_cell,
// or compute the wilber-interp reference partition and score it the same way.
//
//   score_sse bounds SORTED.f64 BOUNDS.txt
//       BOUNDS.txt: k+1 integer bin edges, 0-based, first 0, last n, nondecreasing.
//       prints: SSE=<%.17g> n=<n> k=<k>
//
//   score_sse wilber SORTED.f64 K
//       Runs unweighted Grønlund wilber-interp (interpolation search, same as
//       v10's patched bench/weighted_wilber) once on the already-sorted array
//       and scores the path with raw_sse. prints: SSE=... and a CELL line.
#include "kmeans.h"

#include <sys/resource.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static double peak_rss_mb() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__)
    return static_cast<double>(ru.ru_maxrss) / (1024.0 * 1024.0);  // bytes
#else
    return static_cast<double>(ru.ru_maxrss) / 1024.0;             // KiB
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

// Identical to v10 bench/scale_cell.cpp raw_sse.
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

static std::vector<int> load_bounds(const char* path, int n) {
    FILE* fp = std::fopen(path, "r");
    if (!fp) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(2); }
    std::vector<int> b;
    int v;
    while (std::fscanf(fp, "%d", &v) == 1) b.push_back(v);
    std::fclose(fp);
    if (b.size() < 2) { std::fprintf(stderr, "bounds file %s has fewer than 2 edges\n", path); std::exit(2); }
    if (b.front() != 0 || b.back() != n) {
        std::fprintf(stderr, "bounds must start at 0 and end at n=%d (got %d .. %d)\n",
                     n, b.front(), b.back());
        std::exit(2);
    }
    for (std::size_t i = 1; i < b.size(); ++i) {
        if (b[i] < b[i - 1]) {
            std::fprintf(stderr, "bounds are not nondecreasing at %zu: %d then %d\n", i, b[i - 1], b[i]);
            std::exit(2);
        }
    }
    return b;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
            "usage: score_sse bounds SORTED.f64 BOUNDS.txt\n"
            "       score_sse wilber SORTED.f64 K\n");
        return 2;
    }
    const std::string mode = argv[1];
    const char* path = argv[2];
    std::vector<double> sorted = load_f64(path);
    const std::size_t n = sorted.size();
    // Caller hands a file already sorted by sort_f64 (std::sort, -0 canonical).
    for (double& x : sorted) x = x + 0.0;
    if (!std::is_sorted(sorted.begin(), sorted.end())) {
        std::fprintf(stderr, "input is not sorted: %s\n", path);
        return 2;
    }
    int D = n ? 1 : 0;
    for (std::size_t i = 1; i < n; ++i) if (sorted[i] != sorted[i - 1]) ++D;

    if (mode == "bounds") {
        if (argc < 4) { std::fprintf(stderr, "usage: score_sse bounds SORTED.f64 BOUNDS.txt\n"); return 2; }
        const std::vector<int> b = load_bounds(argv[3], static_cast<int>(n));
        const double sse = raw_sse(sorted, b);
        const int k = static_cast<int>(b.size()) - 1;
        std::printf("SSE=%.17g n=%zu k=%d D=%d\n", sse, n, k, D);
        std::fflush(stdout);
        return 0;
    }

    if (mode == "wilber") {
        if (argc < 4) { std::fprintf(stderr, "usage: score_sse wilber SORTED.f64 K\n"); return 2; }
        const int k = std::atoi(argv[3]);
        if (k < 1) { std::fprintf(stderr, "bad k\n"); return 2; }
        if (static_cast<std::size_t>(k) >= n || k > D) {
            std::printf("CELL method=wilber k=%d n=%zu D=%d status=skip_k_ge_D peak_rss_mb=%.1f\n",
                        k, n, D, peak_rss_mb());
            return 0;
        }
        kmeans_wilber fn(sorted);
        fn.set_search_strategy(search_strategy::INTERPOLATION);
        auto res = fn.compute_and_report(static_cast<std::size_t>(k));
        if (!res) {
            std::fprintf(stderr, "wilber returned null\n");
            return 4;
        }
        std::vector<int> b;
        for (auto it = res->path.rbegin(); it != res->path.rend(); ++it)
            b.push_back(static_cast<int>(*it));
        const double sse = raw_sse(sorted, b);
        const int k_eff = static_cast<int>(res->path.size()) - 1;
        std::printf("SSE=%.17g n=%zu k=%d D=%d k_eff=%d lib_cost=%.17g peak_rss_mb=%.1f\n",
                    sse, n, k, D, k_eff, res->cost, peak_rss_mb());
        std::printf("CELL method=wilber k=%d n=%zu D=%d status=ok sse=%.17g lib_cost=%.17g "
                    "k_eff=%d peak_rss_mb=%.1f\n",
                    k, n, D, sse, res->cost, k_eff, peak_rss_mb());
        std::fflush(stdout);
        return 0;
    }

    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
