// LEX quickstart in C++: three contrasting distributions side by side. Generates its own data.
// For uniform, normal and lognormal samples (N = 200,000, 10 % of values rounded to 2 decimals so
// there are ties) it fits k = 16 bins with lex_hybrid (the default), lex_plain and ckmeans_exact, and prints the
// gap to the exact optimum, the uniform min-max MSE ratio and the time. The random numbers come
// from splitmix64 + Box-Muller (no std:: distributions, whose output differs between standard
// libraries), so the data match the Scala example up to libm rounding of log, cos and exp.
#include <lex/lex.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace {

struct Rng { // splitmix64
    std::uint64_t s;
    std::uint64_t next() {
        std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; } // [0, 1)
    double normal() {
        const double u = 1.0 - uniform(), v = uniform();
        return std::sqrt(-2.0 * std::log(u)) * std::cos(6.283185307179586 * v);
    }
};

std::vector<double> sample(const std::string& dist, Rng& rng, std::size_t n) {
    std::vector<double> x(n);
    for (double& v : x) {
        v = dist == "uniform" ? rng.uniform() : dist == "normal" ? rng.normal() : std::exp(rng.normal());
        if (rng.uniform() < 0.10) v = std::round(v * 100) / 100; // some exact duplicates
    }
    return x;
}

double uniform_minmax_mse(const std::vector<double>& x, int k) {
    const auto [lo, hi] = std::minmax_element(x.begin(), x.end());
    const double step = (*hi - *lo) / k;
    double s = 0;
    for (double v : x) {
        const double cell = std::min(std::floor((v - *lo) / step), k - 1.0);
        const double q = *lo + (cell + 0.5) * step;
        s += (v - q) * (v - q);
    }
    return s / static_cast<double>(x.size());
}

template <class F> double ms(F&& f) {
    const auto t = std::chrono::steady_clock::now();
    f();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}

} // namespace

int main() {
    const std::size_t n = 200000;
    const int k = 16;
    Rng rng{2026};
    std::printf("%-10s %8s %3s %12s %10s %12s %10s %8s %8s\n", "dist", "N", "k", "gap hybrid", "gap plain",
                "uniform/LEX", "t hybrid", "t plain", "t exact");
    for (const std::string dist : {"uniform", "normal", "lognormal"}) {
        const std::vector<double> x = sample(dist, rng, n);
        lex::Result d, p, e;
        const double td = ms([&] { d = lex::lex_hybrid(x, k); });
        const double tp = ms([&] { p = lex::lex_plain(x, k); });
        const double te = ms([&] { e = lex::ckmeans_exact(x, k); });
        std::printf("%-10s %8zu %3d %11.2e%% %9.3f%% %11.2fx %8.1fms %6.1fms %6.1fms\n", dist.c_str(), n, d.k(),
                    100 * (d.sse / e.sse - 1), 100 * (p.sse / e.sse - 1),
                    uniform_minmax_mse(x, k) / (d.sse / static_cast<double>(n)), td, tp, te);
    }
    // weighted: same API with a weight per value
    const std::vector<double> x = sample("lognormal", rng, n);
    std::vector<double> w(n);
    for (double& v : w) v = std::exp(rng.normal());
    const lex::Result rw = lex::lex_hybrid(x, w, k), ew = lex::ckmeans_exact(x, w, k);
    std::printf("weighted lognormal: gap hybrid %.2e%%, k used %d\n", 100 * (rw.sse / ew.sse - 1), rw.k());
    return 0;
}
