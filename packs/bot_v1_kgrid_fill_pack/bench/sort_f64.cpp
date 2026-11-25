// Sort a little-endian float64 file the same way as v10 scale_cell:
//   x = x + 0.0  (canonicalise -0), then std::sort.
// Usage: sort_f64 IN.f64 OUT.f64 DATASET REPS
// One unrecorded warm-up + REPS timed sorts; sort_s is the median.
// Prints one CELL line (sort_s, n, D, peak_rss_mb).
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

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: sort_f64 IN.f64 OUT.f64 DATASET REPS\n");
        return 2;
    }
    const char* path = argv[1];
    const char* outp = argv[2];
    const char* dataset = argv[3];
    const int reps = std::max(1, std::atoi(argv[4]));
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
