// Minimal test support: CHECK macros, a data loader and the O(k d^2) DP oracle.
#pragma once

#include "lex/lex.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <vector>

inline int g_failures = 0;
inline int g_checks = 0;

#define CHECK(cond)                                                                      \
    do {                                                                                 \
        ++g_checks;                                                                      \
        if (!(cond)) { ++g_failures; std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } \
    } while (0)
#define CHECK_MSG(cond, msg)                                                             \
    do {                                                                                 \
        ++g_checks;                                                                      \
        if (!(cond)) { ++g_failures; std::printf("FAIL %s:%d: %s [%s]\n", __FILE__, __LINE__, #cond, std::string(msg).c_str()); } \
    } while (0)
#define CHECK_NEAR(a, b, tol) CHECK_MSG(std::abs((a) - (b)) <= (tol), std::to_string(a) + " vs " + std::to_string(b))

inline void test(const char* name, const std::function<void()>& body) {
    const int before = g_failures;
    body();
    std::printf("%s %s\n", g_failures == before ? "  +" : "  X", name);
}

inline int finish() {
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}

/// The tests name each LEX method by its Options: coarse_groups = 0 runs lex_plain, anything else
/// lex_hybrid. test_coarse checks that this matches calling lex_plain directly.
inline lex::Options plain() { lex::Options o; o.coarse_groups = 0; return o; }
inline lex::Options plain_deck() { lex::Options o = lex::deck(); o.coarse_groups = 0; return o; } // lex_plain(xs, k, deck())
inline lex::Result solve(const std::vector<double>& xs, int k, const lex::Options& o = {}) {
    return o.coarse_groups == 0 ? lex::lex_plain(xs, k, o) : lex::lex_hybrid(xs, k, o);
}
inline lex::Result solve(const std::vector<double>& xs, const std::vector<double>& ws, int k, const lex::Options& o = {}) {
    return o.coarse_groups == 0 ? lex::lex_plain(xs, ws, k, o) : lex::lex_hybrid(xs, ws, k, o);
}

/// One-column file; skips a header line if it does not parse as a number.
inline std::vector<double> load_column(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::printf("missing %s\n", path.c_str()); std::exit(2); }
    std::vector<double> xs;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        char* end = nullptr;
        const double v = std::strtod(line.c_str(), &end);
        if (end == line.c_str()) continue; // header
        xs.push_back(v);
    }
    return xs;
}

/// Number of set bits (portable; __builtin_popcount is GCC/Clang only).
inline int popcount(unsigned v) {
    int c = 0;
    for (; v; v &= v - 1) ++c;
    return c;
}

inline std::vector<double> take(const std::vector<double>& xs, std::size_t n) {
    return {xs.begin(), xs.begin() + static_cast<long>(std::min(n, xs.size()))};
}

struct Optimum { double sse; std::vector<int> counts; };

/// Exact optimum by O(k d^2) DP over distinct values; bin SSEs by weighted Welford updates.
inline Optimum oracle(std::vector<double> s, int k) {
    std::sort(s.begin(), s.end());
    std::vector<double> vals;
    std::vector<int> cnt;
    for (double v : s) {
        if (vals.empty() || v != vals.back()) { vals.push_back(v); cnt.push_back(0); }
        ++cnt.back();
    }
    const int d = static_cast<int>(vals.size());
    std::vector<std::vector<double>> cost(d, std::vector<double>(d + 1));
    for (int i = 0; i < d; ++i) {
        double n = 0, mean = 0, m2 = 0;
        for (int j = i; j < d; ++j) {
            const double c = cnt[j], delta = vals[j] - mean, n2 = n + c;
            mean += delta * c / n2;
            m2 += delta * delta * n * c / n2;
            n = n2;
            cost[i][j + 1] = m2;
        }
    }
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> best(k + 1, std::vector<double>(d + 1, inf));
    std::vector<std::vector<int>> from(k + 1, std::vector<int>(d + 1));
    best[0][0] = 0;
    for (int b = 1; b <= k; ++b)
        for (int j = b; j <= d; ++j)
            for (int i = b - 1; i < j; ++i)
                if (const double c = best[b - 1][i] + cost[i][j]; c < best[b][j]) { best[b][j] = c; from[b][j] = i; }
    std::vector<int> counts(k);
    for (int b = k, j = d; b >= 1; --b) {
        const int i = from[b][j];
        int c = 0;
        for (int t = i; t < j; ++t) c += cnt[t];
        counts[b - 1] = c;
        j = i;
    }
    return {best[k][d], counts};
}

/// Exact SSE of given counts over the sorted values (two passes per bin).
inline double sse_of(std::vector<double> s, const std::vector<int>& counts) {
    std::sort(s.begin(), s.end());
    double total = 0;
    std::size_t start = 0;
    for (int c : counts) {
        double m = 0;
        for (int i = 0; i < c; ++i) m += s[start + i];
        m /= c;
        for (int i = 0; i < c; ++i) total += (s[start + i] - m) * (s[start + i] - m);
        start += c;
    }
    return total;
}

inline int distinct(std::vector<double> xs) {
    std::sort(xs.begin(), xs.end());
    return static_cast<int>(std::unique(xs.begin(), xs.end()) - xs.begin());
}
