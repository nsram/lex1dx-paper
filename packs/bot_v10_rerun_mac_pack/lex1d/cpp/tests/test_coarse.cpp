// Gates for the coarse exact start (Options::coarse_groups):
//  (b) coarse_groups >= atoms starts from the global optimum, so the result matches ckmeans_exact;
//  (c) the DP on M groups never does worse than Collapse to k, and refine never makes it worse;
//  (d) edge cases and option interactions; (e) on tiny inputs the DP on groups is optimal among
//  partitions that never split a group (brute force). Gate (a), coarse_groups = 0 unchanged, is
//  the Scala cross-check plus the "M <= k" case below.
#include "support.hpp"
#include "../src/detail.hpp"

#include <random>
#include <stdexcept>

using namespace lex;
using detail::Atoms;

static Options coarse(std::size_t m, Options o = {}) { o.coarse_groups = m; return o; }

static Atoms atoms_of(std::vector<double> xs, double precision = 0) {
    std::sort(xs.begin(), xs.end());
    return Atoms(std::move(xs), precision);
}

/// Collapse to m groups, then the DP on the groups: the coarse exact start before refinement.
static std::vector<int> dp_on_groups(const Atoms& a, int m, int k, std::vector<int>& groups) {
    groups = detail::collapse(a, m);
    std::vector<double> means(m);
    for (int g = 0; g < m; ++g) means[g] = a.mean(groups[g], groups[g + 1]);
    return detail::exact(groups, means, k);
}

static bool valid(const std::vector<int>& b, int n, int k) {
    if (static_cast<int>(b.size()) != k + 1 || b.front() != 0 || b.back() != n) return false;
    for (int i = 1; i <= k; ++i) if (b[i] <= b[i - 1]) return false;
    return true;
}

static bool same(const Result& x, const Result& y) {
    return x.bounds == y.bounds && x.sse == y.sse && x.centers == y.centers && x.start == y.start &&
           x.sweeps == y.sweeps && x.moves == y.moves && x.atoms == y.atoms;
}

static std::vector<double> random_data(std::mt19937_64& rng, int n, int kind) {
    std::vector<double> xs(n);
    std::normal_distribution<double> g;
    std::lognormal_distribution<double> ln(0, 1);
    for (double& x : xs) {
        switch (kind) {
            case 0: x = g(rng); break;
            case 1: x = ln(rng); break;
            case 2: x = static_cast<double>(rng() % 50); break; // many duplicates
            default: x = (rng() % 2 ? 10 : 0) + g(rng); break;   // bimodal
        }
    }
    return xs;
}

int main() {
    test("coarse_groups >= atoms gives ckmeans_exact's partition and SSE", [] {
        std::mt19937_64 rng(11);
        for (int t = 0; t < 400; ++t) {
            const auto xs = random_data(rng, 2 + static_cast<int>(rng() % 300), t % 4);
            const int d = distinct(xs);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(d, 20)));
            const Result e = ckmeans_exact(xs, k);
            for (std::size_t m : {static_cast<std::size_t>(d), static_cast<std::size_t>(d) + 7, SIZE_MAX}) {
                const Result r = solve(xs, k, coarse(m));
                CHECK(r.sse <= e.sse);
                CHECK_MSG(r.bounds == e.bounds, "t=" + std::to_string(t) + " k=" + std::to_string(k));
                CHECK(r.start == (k < d ? Start::CoarseExact : Start::Collapse));
            }
        }
    });

    test("DP on M groups <= Collapse to k (nested merges); refine never raises it", [] {
        std::mt19937_64 rng(12);
        for (int t = 0; t < 1500; ++t) {
            const double precision = t % 5 == 4 ? 1e-2 : 0;
            const auto xs = random_data(rng, 5 + static_cast<int>(rng() % 2000), t % 4);
            const Atoms a = atoms_of(xs, precision);
            if (a.d < 3) continue;
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(a.d - 2, 64)));
            const int m = k + 1 + static_cast<int>(rng() % static_cast<unsigned>(a.d - k - 1)); // k < m < d
            std::vector<int> groups;
            const std::vector<int> b = dp_on_groups(a, m, k, groups);
            const std::vector<int> c = detail::collapse(a, k);
            CHECK(valid(groups, a.n, m));
            CHECK(valid(b, a.n, k));
            for (int x : c) CHECK(std::binary_search(groups.begin(), groups.end(), x)); // Collapse(k) coarsens Collapse(m)
            for (int x : b) CHECK(std::binary_search(groups.begin(), groups.end(), x)); // DP cuts only between groups
            const double sc = a.sse(c), sb = a.sse(b);
            CHECK_MSG(sb <= sc * (1 + 1e-12) + 1e-300, std::to_string(sb) + " vs " + std::to_string(sc));
            Options o = coarse(static_cast<std::size_t>(m));
            o.precision = precision;
            const Result r = solve(xs, k, o);
            CHECK(r.start == Start::CoarseExact);
            CHECK(r.sse <= sse_of(xs, [&] { std::vector<int> s(k); for (int i = 0; i < k; ++i) s[i] = b[i + 1] - b[i]; return s; }()) * (1 + 1e-12) + 1e-300);
            CHECK(r.sse >= ckmeans_exact(xs, k).sse * (1 - 1e-12) || precision > 0);
        }
    });

    test("DP on groups is optimal among partitions that never split a group (brute force)", [] {
        std::mt19937_64 rng(13);
        for (int t = 0; t < 600; ++t) {
            const int n = 3 + static_cast<int>(rng() % 14);
            std::vector<double> xs(n);
            for (double& x : xs) x = static_cast<double>(rng() % (t % 2 ? 9 : 1000)) / 7.0;
            std::sort(xs.begin(), xs.end());
            const Atoms a = atoms_of(xs);
            if (a.d < 3) continue;
            const int m = 2 + static_cast<int>(rng() % static_cast<unsigned>(a.d - 1)); // 2..d
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(m));      // 1..m
            std::vector<int> groups;
            const std::vector<int> b = dp_on_groups(a, m, k, groups);
            double best = std::numeric_limits<double>::infinity();
            for (unsigned mask = 0; mask < (1u << (m - 1)); ++mask) { // cut after group i if bit i set
                if (popcount(mask) != k - 1) continue;
                std::vector<int> counts;
                int prev = 0;
                for (int i = 0; i < m - 1; ++i)
                    if (mask >> i & 1) { counts.push_back(groups[i + 1] - prev); prev = groups[i + 1]; }
                counts.push_back(n - prev);
                best = std::min(best, sse_of(xs, counts));
            }
            std::vector<int> counts(k);
            for (int i = 0; i < k; ++i) counts[i] = b[i + 1] - b[i];
            CHECK_NEAR(sse_of(xs, counts), best, 1e-9 * (1 + best));
        }
    });

    test("M <= k (coarse_groups <= k) is the plain Collapse start, bit for bit", [] {
        std::mt19937_64 rng(14);
        for (int t = 0; t < 200; ++t) {
            const auto xs = random_data(rng, 10 + static_cast<int>(rng() % 500), t % 4);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(distinct(xs), 16)));
            for (const Options& base : {plain(), plain_deck()}) {
                const Result ref = solve(xs, k, base);
                for (std::size_t m : {std::size_t{1}, static_cast<std::size_t>(k)})
                    CHECK(same(solve(xs, k, coarse(m, base)), ref));
            }
        }
    });

    test("edge cases: n = 0/1, k = 1, k = atoms, all equal, non-finite input", [] {
        const double inf = std::numeric_limits<double>::infinity();
        auto rejects = [](const std::vector<double>& xs, int k, const std::string& fragment) {
            try { solve(xs, k, coarse(64)); CHECK_MSG(false, "no exception for " + fragment); }
            catch (const std::invalid_argument& e) { CHECK_MSG(std::string(e.what()).find(fragment) != std::string::npos, e.what()); }
        };
        rejects({}, 1, "empty");
        rejects({1.0, std::nan("")}, 1, "index 1 is NaN");
        rejects({0.0, inf}, 1, "index 1 is Infinity");
        rejects({1.0, 2.0}, 0, "k must be at least 1");
        CHECK(solve(std::vector<double>{1.0, 1.0, 2.0}, 3, coarse(64)).k() == 2); // clamped to the atoms

        const Result one = solve(std::vector<double>{42.0}, 1, coarse(8));
        CHECK(one.sse == 0.0); CHECK(one.bounds == (std::vector<int>{0, 1})); CHECK(one.start == Start::Collapse);

        const std::vector<double> xs{1, 2, 4, 8, 8, 9, 30};
        const Result k1 = solve(xs, 1, coarse(4));
        CHECK(k1.bounds == (std::vector<int>{0, 7})); CHECK(k1.start == Start::CoarseExact);
        CHECK_NEAR(k1.sse, sse_of(xs, {7}), 1e-12);

        const Result kd = solve(xs, 6, coarse(100)); // k = atoms: nothing to merge
        CHECK(kd.sse == 0.0); CHECK(same(kd, solve(xs, 6)));

        const Result eq = solve(std::vector<double>(9, -2.5), 1, coarse(16));
        CHECK(eq.sse == 0.0); CHECK(eq.atoms == 1); CHECK(eq.bounds == (std::vector<int>{0, 9}));
    });

    test("multi_start with coarse_groups: coarse candidate replaces Collapse; bounds valid", [] {
        std::mt19937_64 rng(15);
        for (int t = 0; t < 200; ++t) {
            const auto xs = random_data(rng, 50 + static_cast<int>(rng() % 3000), t % 4);
            const Atoms a = atoms_of(xs, 1e-4);
            if (a.d < 4) continue;
            const int k = 2 + static_cast<int>(rng() % static_cast<unsigned>(std::min(a.d - 2, 32)));
            const Result r = solve(xs, k, coarse(static_cast<std::size_t>(8 * k), plain_deck()));
            CHECK(r.start != Start::Collapse);
            CHECK(valid(r.bounds, a.n, k));
        }
    });

    test("atoms <= M: lex_hybrid returns ckmeans_exact's partition with 0 refine moves; atoms = M + 1 is not exact", [] {
        // Each dataset with its distinct-value count d. Run with M = d (the D == M edge), M = d + 1 and
        // the default M (4096, >= d except for the last dataset): the DP on all atoms is the optimum,
        // so Refine has nothing to move. With M = d - 1 one Collapse merge happens, and only
        // SSE >= exact is guaranteed.
        std::vector<std::pair<const char*, std::vector<double>>> sets;
        sets.push_back({"evenly spaced ties", {1, 2, 3, 4, 5, 6, 7, 8}}); // many equal-SSE partitions
        sets.push_back({"duplicates, d = 5", {0, 0, 0, 1, 1, 2, 2, 2, 2, 5, 9, 9}});
        std::vector<double> dup(3000), gauss(300), big(4096);
        std::mt19937_64 rng(4096);
        std::normal_distribution<double> g;
        for (double& x : dup) x = static_cast<double>(rng() % 40); // d = 40, ~75 copies each
        for (double& x : gauss) x = g(rng);
        for (std::size_t i = 0; i < big.size(); ++i) big[i] = g(rng);
        sets.push_back({"integers 0..39, n = 3000", dup});
        sets.push_back({"gaussian, d = 300", gauss});
        std::vector<double> big2 = big;                       // d = 4096 = default M, with 4096 duplicates added
        big2.insert(big2.end(), big.begin(), big.end());
        sets.push_back({"gaussian, d = 4096 (default M)", big2});
        for (const auto& [name, xs] : sets) {
            const int d = distinct(xs);
            for (int k : {2, 3, 7, 16, 64}) {
                if (k >= d) continue;
                const Result e = ckmeans_exact(xs, k);
                for (std::size_t m : {static_cast<std::size_t>(d), static_cast<std::size_t>(d) + 1, Options{}.coarse_groups}) {
                    if (m < static_cast<std::size_t>(d)) continue;
                    const Result r = solve(xs, k, coarse(m));
                    const std::string at = std::string(name) + " k=" + std::to_string(k) + " M=" + std::to_string(m);
                    CHECK_MSG(r.start == Start::CoarseExact, at);
                    CHECK_MSG(r.bounds == e.bounds, at);
                    CHECK_MSG(r.sse == e.sse, at);
                    CHECK_MSG(r.moves == 0, at + " moves=" + std::to_string(r.moves));
                }
                const Result c = solve(xs, k, coarse(static_cast<std::size_t>(d) - 1)); // atoms = M + 1
                CHECK_MSG(c.sse >= e.sse * (1 - 1e-12), std::string(name) + " k=" + std::to_string(k));
            }
        }
        std::vector<double> big3 = big; // d = 4097: one atom above the default M
        big3.push_back(10.0);
        for (int k : {2, 16, 64}) CHECK(lex_hybrid(big3, k).sse >= ckmeans_exact(big3, k).sse * (1 - 1e-12));
    });

    test("defaults: lex_hybrid uses coarse 4096; lex_plain ignores coarse_groups; deck()", [] {
        CHECK(Options{}.coarse_groups == 4096);
        CHECK(deck().coarse_groups == 4096 && deck().precision == 1e-4 && deck().multi_start);
        std::mt19937_64 rng(16384);
        std::normal_distribution<double> g;
        std::vector<double> xs(3000), ws(3000);
        for (double& x : xs) x = g(rng);
        for (double& w : ws) w = 0.5 + std::abs(g(rng));
        for (int k : {2, 8, 32}) {
            const Result d = lex_hybrid(xs, k), e = ckmeans_exact(xs, k), p = lex_plain(xs, k);
            CHECK(d.start == Start::CoarseExact);
            CHECK(d.bounds == e.bounds); // atoms <= 4096: the coarse start is the DP itself
            CHECK(p.start == Start::Collapse);
            CHECK(p.sse >= e.sse);
            // lex_plain is lex_hybrid with coarse_groups = 0, whatever coarse_groups it is given
            CHECK(same(p, lex_hybrid(xs, k, plain())));
            CHECK(same(lex_plain(xs, k, deck()), lex_hybrid(xs, k, plain_deck())));
            CHECK(same(lex_plain(xs, ws, k, deck()), lex_hybrid(xs, ws, k, plain_deck())));
            CHECK(same(lex_plain(xs, ws, k), lex_hybrid(xs, ws, k, plain())));
        }
    });

    return finish();
}
