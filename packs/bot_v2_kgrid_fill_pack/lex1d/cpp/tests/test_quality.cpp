#include "support.hpp"
#include "../src/detail.hpp"

#include <random>

using namespace lex;

int main() {
    test("random small instances vs exact DP optimum (report; LEX is a local optimum)", [] {
        std::mt19937_64 rng(20261112);
        auto uni = [&](int m) { return static_cast<int>(rng() % static_cast<unsigned>(m)); };
        for (bool multi : {false, true}) {
            int cases = 0, optimal = 0;
            double worst = 0;
            for (int t = 0; t < 3000; ++t) {
                std::vector<double> xs(3 + uni(28));
                for (double& x : xs) x = uni(21);
                const int k = 1 + uni(std::min(6, distinct(xs)));
                Options o = plain(); o.multi_start = multi;
                const Result r = solve(xs, k, o);
                const double opt = oracle(xs, k).sse;
                const double sse = sse_of(xs, r.sizes);
                CHECK_NEAR(r.sse, sse, 1e-9 * (1 + sse));
                CHECK(sse >= opt - 1e-9 * (1 + opt));
                ++cases;
                if (sse <= opt + 1e-9 * (1 + opt)) ++optimal;
                else worst = std::max(worst, (sse - opt) / opt);
            }
            std::printf("ORACLE multi_start=%d: %d/%d optimal (%.1f%%), worst relative gap %.1f%%\n", multi, optimal, cases,
                        100.0 * optimal / cases, 100 * worst);
        }
    });

    test("every applied move lowers the exact SSE (tie moves keep it); atoms are never split", [] {
        std::mt19937_64 rng(7);
        auto uni = [&](int m) { return static_cast<int>(rng() % static_cast<unsigned>(m)); };
        auto check = [](std::vector<double> xs, int k, bool multi) {
            std::sort(xs.begin(), xs.end());
            const detail::Atoms a(xs, 0);
            std::vector<int> b = multi ? detail::equal_frequency(a.n, k) : detail::collapse(a, k);
            double last = a.sse(b);
            bool ok = true;
            std::vector<int> prev = b;
            detail::refine(a, b, 1e-11, [&](const std::vector<int>& bb) {
                // gain moves strictly lower SSE; tie moves (bounds move right) keep it up to rounding
                const double s = a.sse(bb);
                const bool right = std::lexicographical_compare(prev.rbegin(), prev.rend(), bb.rbegin(), bb.rend());
                if (!(s < last || (right && s <= last + 1e-12 * (1 + std::abs(last))))) ok = false;
                last = s;
                prev = bb;
            });
            CHECK(ok);
            if (!multi)
                for (int e : b) CHECK(e == a.n || a.start[a.atom_id[e]] == e);
        };
        for (int t = 0; t < 500; ++t) {
            std::vector<double> xs(5 + uni(60));
            for (double& x : xs) x = uni(40);
            const int d = distinct(xs);
            check(xs, 1 + uni(std::min(8, d)), false);
            check(xs, 1 + uni(std::min(8, d)), true);
        }
        std::normal_distribution<double> g;
        std::vector<double> big(20000);
        for (double& x : big) x = g(rng);
        check(big, 16, false);
    });

    test("translation invariance: shifts up to 1e8 give the same SSE, no crash or hang", [] {
        std::mt19937_64 rng(42);
        std::normal_distribution<double> g;
        std::vector<double> base(50000);
        for (double& x : base) x = g(rng);
        const Result ref = solve(base, 16);
        for (double shift : {1e3, 1e6, 1e8}) {
            auto xs = base;
            for (double& x : xs) x += shift;
            CHECK(std::abs(solve(xs, 16).sse - ref.sse) / ref.sse < 1e-6);
        }
        auto xs = base;
        for (double& x : xs) x += 1e6;
        const Options grouped{1e-4};
        const double g0 = solve(base, 16, grouped).sse, g6 = solve(xs, 16, grouped).sse;
        CHECK(std::abs(g6 - g0) / g0 < 1e-6);
    });

    test("2024 deck config reproduces the deck's printed MSE (LEX_SIM; skipped if unset)", [] {
        const char* sim = std::getenv("LEX_SIM");
        if (!sim) { std::printf("    skipped: set LEX_SIM to lex/simulation/input\n"); return; }
        struct Case { const char* dist; double printed; double old; bool same_path; int digits; };
        const Case cases[] = {{"uniform_continuous", 0.00002032, 0.000020320495, true, 4},
                              {"gaussian_continuous", 0.0006447, 0.000644686140, true, 4},
                              {"exponential_continuous", 0.0005397, 0.000539735961, false, 4}};
        for (const Case& c : cases) {
            const auto all = load_column(std::string(sim) + "/" + c.dist + ".txt");
            const std::vector<double> xs(all.end() - 1048576, all.end());
            const Result r = solve(xs, 64, plain_deck());
            const double mse = r.sse / 1048576;
            char a[32], b[32];
            std::snprintf(a, sizeof a, "%.*e", c.digits - 1, mse);
            std::snprintf(b, sizeof b, "%.*e", c.digits - 1, c.printed);
            std::printf("    DECK|%s|mse=%.12f|start=%d|atoms=%d|sweeps=%d|moves=%d\n", c.dist, mse, static_cast<int>(r.start), r.atoms,
                        r.sweeps, r.moves);
            CHECK_MSG(std::string(a) == b, std::string(a) + " vs " + b);
            if (c.same_path) CHECK_NEAR(mse, c.old, 5e-13);
        }
    });

    return finish();
}
