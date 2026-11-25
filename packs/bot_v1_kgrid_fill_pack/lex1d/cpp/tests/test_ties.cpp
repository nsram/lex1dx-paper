// Tie-break rule and the sorted-input fast path.
//  (a) equal values always share an atom (any precision), so no bin bound splits equal values;
//  (b) results are bit-identical under any permutation of the input (also weighted, with equal
//      values carrying different weights), and -0.0 is treated as +0.0;
//  (c) small exact-tie cases: the DP keeps the largest cut index and Refine the move with the
//      largest new bounds, so ties resolve to the largest cuts as in Ckmeans.1d.dp (R);
//  (d) already-sorted input skips the sort (Result::presorted) and gives the same result.
#include "support.hpp"
#include "../src/detail.hpp"

#include <cstring>
#include <random>

using namespace lex;
using detail::Atoms;

static bool same(const Result& x, const Result& y) {
    auto bits = [](const std::vector<double>& a, const std::vector<double>& b) {
        return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
    };
    return x.bounds == y.bounds && x.sizes == y.sizes && std::memcmp(&x.sse, &y.sse, sizeof(double)) == 0 &&
           bits(x.centers, y.centers) && bits(x.upper, y.upper) && x.start == y.start && x.sweeps == y.sweeps &&
           x.moves == y.moves && x.atoms == y.atoms;
}

static Options coarse(std::size_t m, Options o = {}) { o.coarse_groups = m; return o; }

/// Heavy-tie data: kind 0 small integers, 1 Gaussian rounded to 0.1, 2 rounded to 1e-4 (deck's
/// grid), 3 a few values repeated thousands of times plus singletons.
static std::vector<double> tied(std::mt19937_64& rng, int kind, int n) {
    std::vector<double> xs(n);
    std::normal_distribution<double> g(0.0, 1.0);
    std::uniform_int_distribution<int> u(0, 20);
    for (int i = 0; i < n; ++i) {
        switch (kind) {
            case 0: xs[i] = u(rng); break;
            case 1: xs[i] = std::round(g(rng) * 10) / 10; break;
            case 2: xs[i] = std::round(g(rng) * 1e4) / 1e4; break;
            default: xs[i] = (i % 7 == 0) ? g(rng) : static_cast<double>(u(rng) % 4); break;
        }
    }
    return xs;
}

static void sorted_of(std::vector<double>& s) { std::sort(s.begin(), s.end()); }

int main() {
    std::mt19937_64 rng(20261114);

    test("equal values share an atom at every precision", [&] {
        for (int kind = 0; kind < 4; ++kind)
            for (double p : {0.0, 1e-4, 0.05, 0.3, 2.5}) {
                std::vector<double> s = tied(rng, kind, 3000);
                sorted_of(s);
                const Atoms a(s, p);
                bool ok = true;
                for (int j = 1; j < a.d; ++j) ok = ok && a.sorted[a.start[j]] != a.sorted[a.start[j] - 1];
                CHECK_MSG(ok, "kind " + std::to_string(kind) + " precision " + std::to_string(p));
            }
        // a run whose mean lags far behind: the precision test alone would still keep equal values
        // together here, and the explicit equality makes it a guarantee
        std::vector<double> s(200, 0.0);
        for (int i = 0; i < 200; ++i) s.push_back(0.3);
        const Atoms a(s, 0.3);
        CHECK(a.d == 2 && a.start[1] == 200);
    });

    test("no bin bound splits equal values", [&] {
        for (int kind = 0; kind < 4; ++kind) {
            const std::vector<double> xs = tied(rng, kind, 5000);
            std::vector<double> s = xs;
            sorted_of(s);
            for (int k : {2, 3, 5, 8}) {
                if (k > distinct(xs)) continue;
                for (const Options& o : {Options{}, plain(), plain_deck(), coarse(6), coarse(1000), coarse(6, plain_deck())}) {
                    const Result r = solve(xs, k, o); // k above the groups is clamped, never an error
                    bool ok = true;
                    for (int i = 1; i < k; ++i) ok = ok && s[r.bounds[i]] != s[r.bounds[i] - 1];
                    CHECK_MSG(ok, "fit kind " + std::to_string(kind) + " k " + std::to_string(k));
                }
                const Result e = ckmeans_exact(xs, k);
                bool ok = true;
                for (int i = 1; i < k; ++i) ok = ok && s[e.bounds[i]] != s[e.bounds[i] - 1];
                CHECK_MSG(ok, "exact kind " + std::to_string(kind) + " k " + std::to_string(k));
            }
        }
    });

    test("results are bit-identical under permutation (unweighted)", [&] {
        for (int kind = 0; kind < 4; ++kind) {
            std::vector<double> xs = tied(rng, kind, 4000);
            const int d = distinct(xs);
            for (int k : {3, 8, 16}) {
                if (k > d) continue;
                std::vector<std::vector<double>> perms{xs};
                std::vector<double> s = xs; sorted_of(s); perms.push_back(s);
                std::reverse(s.begin(), s.end()); perms.push_back(s);
                for (int t = 0; t < 3; ++t) { std::shuffle(xs.begin(), xs.end(), rng); perms.push_back(xs); }
                for (const Options& o : {Options{}, plain(), plain_deck(), coarse(10), coarse(10, plain_deck())}) {
                    const Result r0 = solve(perms[0], k, o);
                    bool ok = true;
                    for (const auto& p : perms) ok = ok && same(solve(p, k, o), r0);
                    CHECK_MSG(ok, "fit kind " + std::to_string(kind) + " k " + std::to_string(k));
                }
                const Result e0 = ckmeans_exact(perms[0], k);
                bool ok = true;
                for (const auto& p : perms) ok = ok && same(ckmeans_exact(p, k), e0);
                CHECK_MSG(ok, "exact kind " + std::to_string(kind) + " k " + std::to_string(k));
            }
        }
    });

    test("results are bit-identical under permutation (weighted, tied values with different weights)", [&] {
        std::uniform_real_distribution<double> uw(0.0, 3.0);
        for (int kind = 0; kind < 4; ++kind) {
            std::vector<double> xs = tied(rng, kind, 3000), ws(xs.size());
            for (double& w : ws) w = uw(rng);
            const int d = distinct(xs);
            std::vector<std::size_t> idx(xs.size());
            for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
            for (int k : {3, 8}) {
                if (k > d) continue;
                for (const Options& o : {Options{}, plain(), plain_deck(), coarse(10)}) {
                    const Result r0 = solve(xs, ws, k, o);
                    bool ok = true;
                    for (int t = 0; t < 4; ++t) {
                        std::shuffle(idx.begin(), idx.end(), rng);
                        std::vector<double> px(xs.size()), pw(xs.size());
                        for (std::size_t i = 0; i < idx.size(); ++i) { px[i] = xs[idx[i]]; pw[i] = ws[idx[i]]; }
                        ok = ok && same(solve(px, pw, k, o), r0);
                        if (o.precision == 0) ok = ok && same(ckmeans_exact(px, pw, k), ckmeans_exact(xs, ws, k));
                    }
                    CHECK_MSG(ok, "weighted kind " + std::to_string(kind) + " k " + std::to_string(k));
                }
            }
        }
    });

    test("-0.0 and +0.0 are one value", [&] {
        std::vector<double> plus, mixed;
        for (int i = 0; i < 300; ++i) {
            const double v = (i % 5) - 2.0; // -2..2, a fifth of them zero
            plus.push_back(v == 0 ? 0.0 : v);
            mixed.push_back(v == 0 ? (i % 2 ? -0.0 : 0.0) : v);
        }
        for (const Options& o : {Options{}, plain(), plain_deck()}) {
            const Result a = solve(plus, 3, o), b = solve(mixed, 3, o);
            CHECK(same(a, b));
            CHECK(a.atoms == 5);
            bool pos = true;
            for (double u : b.upper) pos = pos && !(u == 0 && std::signbit(u));
            CHECK(pos);
        }
        CHECK(same(ckmeans_exact(plus, 3), ckmeans_exact(mixed, 3)));
    });

    test("exact ties: largest cuts, as Ckmeans.1d.dp", [&] {
        // three equally spaced values, k = 2: both partitions have SSE 0.5; R gives sizes 2|1
        CHECK((ckmeans_exact({0, 1, 2}, 2).bounds == std::vector<int>{0, 2, 3})); // DP: largest cut index
        CHECK((solve({0, 1, 2}, 2).bounds == std::vector<int>{0, 2, 3}));
        CHECK((solve({0, 1, 2}, 2, plain()).bounds == std::vector<int>{0, 2, 3}));       // Collapse merges (0,1); no strict gain
        CHECK((ckmeans_exact({2, 2, 0, 1, 1, 0}, 2).bounds == std::vector<int>{0, 4, 6}));
        CHECK((solve({2, 2, 0, 1, 1, 0}, 2).bounds == std::vector<int>{0, 4, 6}));
        // five equally spaced values, k = 4: every adjacent merge costs exactly 0.5
        {
            const Atoms a({0, 1, 2, 3, 4}, 0);
            CHECK((detail::collapse(a, 4) == std::vector<int>{0, 2, 3, 4, 5}));
            CHECK((detail::exact(a, 4) == std::vector<int>{0, 2, 3, 4, 5}));
        }
        // ten values 0..9, k = 3: 4|3|3, 3|4|3 and 3|3|4 all have SSE 9; R gives 4|3|3
        {
            std::vector<double> ten;
            for (int i = 0; i < 10; ++i) ten.push_back(i);
            const std::vector<int> want{0, 4, 7, 10};
            CHECK(ckmeans_exact(ten, 3).bounds == want);
            for (const Options& o : {Options{}, plain(), plain_deck(), coarse(6), coarse(10)}) CHECK(solve(ten, 3, o).bounds == want);
        }
        // Refine from {0},{1,2,3},{4}: moves to {0,1},{2},{3,4}, {0,1},{2,3},{4} and {0},{1,2},{3,4}
        // all reach SSE 1.0; the larger bounds (right bound first) win: {0,1},{2,3},{4}.
        {
            const Atoms a({0, 1, 2, 3, 4}, 0);
            std::vector<int> b{0, 1, 4, 5};
            const auto st = detail::refine(a, b, 1e-11);
            CHECK((b == std::vector<int>{0, 2, 4, 5}));
            CHECK(st.moves == 1);
        }
    });

    test("Refine's unweighted tie window is needed: grid ties end on the largest cuts", [&] {
        // 0..99 ten times, deck preset (lex_plain), k = 8: several partitions have SSE exactly 13000,
        // but the tie move's computed delta is a few ulps off 0, so it is a tie only within
        // kTieEps * mag. With the window set to 0 the result stays on 130|120|130|... cuts.
        std::vector<double> grid;
        for (int i = 0; i < 1000; ++i) grid.push_back(i % 100);
        const std::vector<int> want{0, 130, 260, 390, 520, 640, 760, 880, 1000};
        CHECK(solve(grid, 8, plain_deck()).bounds == want);
        CHECK(ckmeans_exact(grid, 8).bounds == want);
        CHECK(detail::kTieEps > 0);
    });

    test("sorted input skips the sort and gives the same result", [&] {
        std::normal_distribution<double> g(0.0, 1.0);
        std::vector<double> xs(20000);
        for (double& x : xs) x = g(rng);
        std::vector<double> s = xs; sorted_of(s);
        for (const Options& o : {Options{}, plain(), plain_deck(), coarse(2000)}) {
            const Result a = solve(xs, 16, o), b = solve(s, 16, o);
            CHECK(!a.presorted && b.presorted);
            CHECK(same(a, b));
        }
        CHECK(ckmeans_exact(s, 16).presorted && !ckmeans_exact(xs, 16).presorted);
        CHECK(same(ckmeans_exact(s, 16), ckmeans_exact(xs, 16)));
        std::vector<double> desc = s;
        std::reverse(desc.begin(), desc.end());
        CHECK(!solve(desc, 16).presorted);
        CHECK(solve(std::vector<double>(50, 1.5), 1).presorted); // all equal counts as sorted
        // weighted: sorted means sorted by (value, weight)
        std::vector<double> x2{0, 1, 1, 2, 3}, w_up{1, 1, 2, 1, 1}, w_down{1, 2, 1, 1, 1};
        CHECK(solve(x2, w_up, 2).presorted);
        CHECK(!solve(x2, w_down, 2).presorted);
        CHECK(same(solve(x2, w_down, 2), solve(std::vector<double>{1, 0, 1, 2, 3}, std::vector<double>{2, 1, 1, 1, 1}, 2)));
    });

    return finish();
}
