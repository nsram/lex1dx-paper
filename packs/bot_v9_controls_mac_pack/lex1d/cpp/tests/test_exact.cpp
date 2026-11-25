#include "support.hpp"

#include <random>

using namespace lex;

int main() {
    test("ckmeans_exact equals the O(k d^2) DP oracle on 3000 random small instances", [] {
        std::mt19937_64 rng(99);
        auto uni = [&](int m) { return static_cast<int>(rng() % static_cast<unsigned>(m)); };
        for (int t = 0; t < 3000; ++t) {
            std::vector<double> xs(1 + uni(40));
            for (double& x : xs) x = uni(25) - 12.0;
            const int k = 1 + uni(std::min(8, distinct(xs)));
            const Result r = ckmeans_exact(xs, k);
            const Optimum o = oracle(xs, k);
            CHECK_NEAR(r.sse, o.sse, 1e-9 * (1 + o.sse));
            CHECK(r.start == Start::Exact);
        }
    });

    test("ckmeans_exact equals brute force over all partitions on tiny inputs", [] {
        std::mt19937_64 rng(5);
        for (int t = 0; t < 300; ++t) {
            const int n = 1 + static_cast<int>(rng() % 10);
            std::vector<double> xs(n);
            for (double& x : xs) x = static_cast<double>(rng() % 1000) / 7.0;
            std::sort(xs.begin(), xs.end());
            xs.erase(std::unique(xs.begin(), xs.end()), xs.end()); // distinct, so every cut set is valid
            const int m = static_cast<int>(xs.size());
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(m));
            double best = std::numeric_limits<double>::infinity();
            for (unsigned mask = 0; mask < (1u << (m - 1)); ++mask) { // cut after position i if bit i is set
                if (popcount(mask) != k - 1) continue;
                std::vector<int> counts;
                int c = 1;
                for (int i = 0; i < m - 1; ++i, ++c)
                    if (mask >> i & 1) { counts.push_back(c); c = 0; }
                counts.push_back(c);
                best = std::min(best, sse_of(xs, counts));
            }
            CHECK_NEAR(ckmeans_exact(xs, k).sse, best, 1e-9 * (1 + best));
        }
    });

    test("ckmeans_exact matches the oracle on rounded and heavy-tailed generated data, k = 2..12", [] {
        // stand-ins for timing-like data (2,400 values rounded to 0.01, many ties) and a skewed
        // population-like list (800 values, lognormal); generated, so no external data is needed
        std::mt19937_64 rng(2379);
        std::normal_distribution<double> g;
        std::vector<double> times(2400), pops(800);
        for (double& x : times) x = std::round((10.2 + 0.15 * g(rng)) * 100) / 100;
        for (double& x : pops) x = std::round(std::exp(10 + 2 * g(rng)));
        for (const auto* xs : {&times, &pops})
            for (int k = 2; k <= 12; ++k) {
                const Optimum o = oracle(*xs, k);
                const Result r = ckmeans_exact(*xs, k);
                CHECK_NEAR(r.sse, o.sse, 1e-9 * o.sse);
            }
        for (int n : {20, 43}) // short prefixes: SSE equals the oracle optimum to 1e-12
            for (int k = 2; k <= 8; ++k) {
                const std::vector<double> xs(times.begin(), times.begin() + n);
                if (k > distinct(xs)) continue;
                CHECK_NEAR(ckmeans_exact(xs, k).sse, oracle(xs, k).sse, 1e-12);
            }
    });

    test("ckmeans_exact is shift invariant and never above LEX", [] {
        std::mt19937_64 rng(3);
        std::normal_distribution<double> g;
        std::vector<double> base(5000);
        for (double& x : base) x = g(rng);
        const Result ref = ckmeans_exact(base, 16);
        for (double shift : {1e3, 1e6, 1e8}) {
            auto xs = base;
            for (double& x : xs) x += shift;
            const Result r = ckmeans_exact(xs, 16);
            CHECK_MSG(std::abs(r.sse - ref.sse) / ref.sse < 1e-6, "shift " + std::to_string(shift));
        }
        CHECK(solve(base, 16).sse >= ref.sse * (1 - 1e-12));
        CHECK(solve(base, 16, plain_deck()).sse >= ref.sse * (1 - 1e-12));
    });

    return finish();
}
