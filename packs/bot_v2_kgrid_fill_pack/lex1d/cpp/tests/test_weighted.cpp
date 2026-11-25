#include "support.hpp"

#include <random>
#include <stdexcept>

using namespace lex;

/// O(k d^2) weighted DP over distinct values. Zero-weight values are kept as atoms (cost 0).
inline Optimum weighted_oracle(std::vector<double> x, std::vector<double> w, int k) {
    const int n = static_cast<int>(x.size());
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });
    std::vector<double> vals;
    std::vector<double> wt;
    for (int i : idx) {
        if (vals.empty() || x[i] != vals.back()) { vals.push_back(x[i]); wt.push_back(0); }
        wt.back() += w[i];
    }
    const int d = static_cast<int>(vals.size());
    std::vector<std::vector<double>> cost(d, std::vector<double>(d + 1));
    // Two-pass centred bin SSE: slow (O(d^3)) but accurate even with 1e12 vs 1e-12 weights,
    // where the raw-sums form Q - S^2/W loses several digits to cancellation.
    for (int i = 0; i < d; ++i)
        for (int j = i; j < d; ++j) {
            double W = 0, S = 0;
            for (int t = i; t <= j; ++t) { W += wt[t]; S += wt[t] * vals[t]; }
            double c = 0;
            if (W > 0) {
                const double m = S / W;
                for (int t = i; t <= j; ++t) { const double e = vals[t] - m; c += wt[t] * e * e; }
            }
            cost[i][j + 1] = c;
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
        counts[b - 1] = j - i; // atom counts, not used for SSE checks
        j = i;
    }
    return {best[k][d], counts};
}

inline double weighted_sse_of_bounds(std::vector<double> x, std::vector<double> w, const std::vector<int>& bounds) {
    const int n = static_cast<int>(x.size());
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });
    std::vector<double> xs(n), ws(n);
    for (int i = 0; i < n; ++i) { xs[i] = x[idx[i]]; ws[i] = w[idx[i]]; }
    double total = 0;
    for (std::size_t b = 0; b + 1 < bounds.size(); ++b) {
        const int lo = bounds[b], hi = bounds[b + 1];
        double W = 0, S = 0;
        for (int i = lo; i < hi; ++i) { W += ws[i]; S += ws[i] * xs[i]; }
        if (W == 0) continue;
        const double m = S / W;
        for (int i = lo; i < hi; ++i) { const double d = xs[i] - m; total += ws[i] * d * d; }
    }
    return total;
}

int main() {
    test("weighted ckmeans_exact matches the O(k d^2) oracle on random weights", [] {
        std::mt19937_64 rng(20261113);
        std::uniform_real_distribution<double> uni(0.0, 5.0);
        std::normal_distribution<double> g;
        for (int t = 0; t < 400; ++t) {
            const int n = 3 + static_cast<int>(rng() % 30);
            std::vector<double> xs(n), w(n);
            for (int i = 0; i < n; ++i) { xs[i] = g(rng); w[i] = 0.05 + uni(rng); }
            const int d = distinct(xs);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(6, d)));
            const Result r = ckmeans_exact(xs, w, k);
            const Optimum o = weighted_oracle(xs, w, k);
            CHECK_NEAR(r.sse, o.sse, 1e-8 * (1 + o.sse));
            CHECK(r.start == Start::Exact);
            CHECK_NEAR(r.sse, weighted_sse_of_bounds(xs, w, r.bounds), 1e-8 * (1 + r.sse));
        }
    });

    test("weighted ckmeans_exact matches brute force on tiny inputs", [] {
        std::mt19937_64 rng(11);
        std::uniform_real_distribution<double> pos(0.1, 4.0);
        for (int t = 0; t < 80; ++t) {
            const int n = 2 + static_cast<int>(rng() % 8);
            std::vector<double> xs(n), w(n);
            for (int i = 0; i < n; ++i) {
                xs[i] = static_cast<double>(rng() % 20);
                w[i] = pos(rng);
            }
            // unique values so every cut set is valid
            std::sort(xs.begin(), xs.end());
            for (int i = 1; i < n; ++i) if (xs[i] <= xs[i - 1]) xs[i] = xs[i - 1] + 1;
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(n));
            double best = std::numeric_limits<double>::infinity();
            for (unsigned mask = 0; mask < (1u << (n - 1)); ++mask) {
                if (popcount(mask) != k - 1) continue;
                std::vector<int> b(1, 0);
                int c = 1;
                for (int i = 0; i < n - 1; ++i, ++c)
                    if (mask >> i & 1) { b.push_back(b.back() + c); c = 0; }
                b.push_back(n);
                best = std::min(best, weighted_sse_of_bounds(xs, w, b));
            }
            CHECK_NEAR(ckmeans_exact(xs, w, k).sse, best, 1e-8 * (1 + best));
        }
    });

    test("zero weights: exact SSE matches the problem with those points dropped", [] {
        std::mt19937_64 rng(21);
        std::normal_distribution<double> g;
        std::uniform_real_distribution<double> pos(0.2, 3.0);
        for (int t = 0; t < 80; ++t) {
            const int npos = 6 + static_cast<int>(rng() % 12);
            std::vector<double> xs, w;
            for (int i = 0; i < npos; ++i) { xs.push_back(g(rng)); w.push_back(pos(rng)); }
            const int nzero = 1 + static_cast<int>(rng() % 8);
            for (int i = 0; i < nzero; ++i) { xs.push_back(g(rng)); w.push_back(0.0); }
            const int dpos = distinct({xs.begin(), xs.begin() + npos});
            if (dpos < 2) continue;
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(5, dpos)));
            // k must also be <= total distinct
            if (k > distinct(xs)) continue;
            const Result r = ckmeans_exact(xs, w, k);
            const Result s = ckmeans_exact(std::vector<double>(xs.begin(), xs.begin() + npos),
                                       std::vector<double>(w.begin(), w.begin() + npos), k);
            CHECK_NEAR(r.sse, s.sse, 1e-8 * (1 + s.sse));
        }
    });

    test("tiny weights and huge ratios still match the oracle", [] {
        std::mt19937_64 rng(33);
        std::normal_distribution<double> g;
        for (int t = 0; t < 80; ++t) {
            const int n = 8 + static_cast<int>(rng() % 16);
            std::vector<double> xs(n), w(n);
            for (int i = 0; i < n; ++i) {
                xs[i] = g(rng);
                const int kind = static_cast<int>(rng() % 4);
                if (kind == 0) w[i] = 1e-12;
                else if (kind == 1) w[i] = 1e12;
                else if (kind == 2) w[i] = 0.0;
                else w[i] = 1.0;
            }
            if (std::all_of(w.begin(), w.end(), [](double v) { return v == 0; })) w[0] = 1;
            const int d = distinct(xs);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(5, d)));
            const Result r = ckmeans_exact(xs, w, k);
            const double raw = weighted_sse_of_bounds(xs, w, r.bounds);
            CHECK_NEAR(r.sse, raw, 1e-8 * (1 + raw));
            const Optimum o = weighted_oracle(xs, w, k);
            // 1e12 vs 1e-12 is beyond comfortable double dynamic range (~1e16);
            // we still require the same partition quality to about 0.2%.
            const double scale = std::max({1.0, o.sse, r.sse});
            CHECK_MSG(std::abs(r.sse - o.sse) <= 2e-3 * scale,
                      "huge-ratio exact vs oracle " + std::to_string(r.sse) + " vs " + std::to_string(o.sse));
        }
    });

    test("extreme weight ratios terminate (Refine cycle guard)", [] {
        // Found by a 300-seed stress run: before the guard, Refine flipped between two partitions forever.
        const std::vector<double> xs{1.4153503275453148, 0.45275334575251319, 1.6316353555349987, -1.3244473157460113, -0.80571464886831112, 1.4002463070789695, 0.35778628864082129, -0.18037482442377539, -0.24516893024270886, -0.032828437692017064, 0.15315063327204761, 0.37458003340011492, -0.74123225585593044, -0.97659650125445463};
        const std::vector<double> w{1e-12, 1e-12, 1e12, 1e12, 1e-12, 0.0, 1e-12, 0.0, 1e12, 0.0, 0.0, 0.0, 1e-12, 0.0};
        const Result r = ckmeans_exact(xs, w, 4);
        CHECK(r.k() == 4);
        CHECK_NEAR(r.sse, weighted_sse_of_bounds(xs, w, r.bounds), 1e-8 * (1 + r.sse));
        CHECK(solve(xs, w, 4).k() == 4);
    });

    test("moderate weight ratios (1e-6 .. 1e6) match the oracle tightly", [] {
        std::mt19937_64 rng(34);
        std::normal_distribution<double> g;
        std::uniform_real_distribution<double> logw(-6.0, 6.0);
        for (int t = 0; t < 80; ++t) {
            const int n = 8 + static_cast<int>(rng() % 16);
            std::vector<double> xs(n), w(n);
            for (int i = 0; i < n; ++i) { xs[i] = g(rng); w[i] = std::pow(10.0, logw(rng)); }
            const int d = distinct(xs);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(5, d)));
            const Result r = ckmeans_exact(xs, w, k);
            const Optimum o = weighted_oracle(xs, w, k);
            CHECK_NEAR(r.sse, o.sse, 1e-7 * (1 + o.sse));
        }
    });

    test("weighted LEX is never below weighted exact; coarse matches exact on short rows", [] {
        std::mt19937_64 rng(44);
        std::normal_distribution<double> g;
        std::uniform_real_distribution<double> pos(0.1, 5.0);
        for (int t = 0; t < 200; ++t) {
            const int n = 20 + static_cast<int>(rng() % 80);
            std::vector<double> xs(n), w(n);
            for (int i = 0; i < n; ++i) { xs[i] = g(rng); w[i] = pos(rng); }
            const int d = distinct(xs);
            const int k = 1 + static_cast<int>(rng() % static_cast<unsigned>(std::min(8, d)));
            const Result ex = ckmeans_exact(xs, w, k);
            const Result lx = solve(xs, w, k, plain());
            CHECK_MSG(lx.sse >= ex.sse * (1 - 1e-10) - 1e-15, "LEX below exact");
            Options o; o.coarse_groups = 16384;
            const Result c = solve(xs, w, k, o);
            CHECK_NEAR(c.sse, ex.sse, 1e-8 * (1 + ex.sse)); // n is short, so coarse is exact
            Options dopt = plain_deck();
            const Result dk = solve(xs, w, k, dopt);
            // deck groups at 1e-4, so it can be slightly above exact; it must not beat exact on raw SSE
            CHECK_MSG(dk.sse >= ex.sse * (1 - 1e-9) - 1e-15, "deck below exact");
        }
    });

    test("null weights pointer dispatches to the unweighted solver", [] {
        const std::vector<double> xs = {1, 2, 3, 4, 5, 6, 7, 8};
        const Result a = ckmeans_exact(xs, 3);
        const Result b = ckmeans_exact(xs.data(), xs.size(), 3, nullptr);
        CHECK(a.sse == b.sse);
        CHECK(a.bounds == b.bounds);
        const Result c = lex_hybrid(xs, 3);
        const Result d = lex_hybrid(xs.data(), xs.size(), 3, nullptr, Options{});
        CHECK(c.sse == d.sse);
        CHECK(c.bounds == d.bounds);
        const Result e = lex_plain(xs, 3);
        const Result f = lex_plain(xs.data(), xs.size(), 3, nullptr, Options{});
        CHECK(e.sse == f.sse);
        CHECK(e.bounds == f.bounds);
    });

    test("rejects negative, NaN, and all-zero weights; size mismatch throws", [] {
        const std::vector<double> xs = {0, 1, 2, 3};
        try { ckmeans_exact(xs, std::vector<double>{1, -1, 1, 1}, 2); CHECK(false); }
        catch (const std::invalid_argument&) { CHECK(true); }
        try { ckmeans_exact(xs, std::vector<double>{1, 1, NAN, 1}, 2); CHECK(false); }
        catch (const std::invalid_argument&) { CHECK(true); }
        try { ckmeans_exact(xs, std::vector<double>{0, 0, 0, 0}, 2); CHECK(false); }
        catch (const std::invalid_argument&) { CHECK(true); }
        try { ckmeans_exact(xs, std::vector<double>{1, 1}, 2); CHECK(false); }
        catch (const std::invalid_argument&) { CHECK(true); }
    });

    test("unit weights give the same SSE as unweighted, within a tight tolerance", [] {
        std::mt19937_64 rng(7);
        std::normal_distribution<double> g;
        std::vector<double> xs(200);
        for (double& x : xs) x = g(rng);
        std::vector<double> w(xs.size(), 1.0);
        for (int k : {2, 5, 12}) {
            CHECK_NEAR(ckmeans_exact(xs, w, k).sse, ckmeans_exact(xs, k).sse, 1e-9 * (1 + ckmeans_exact(xs, k).sse));
            CHECK_NEAR(solve(xs, w, k).sse, solve(xs, k).sse, 1e-9 * (1 + solve(xs, k).sse));
        }
    });

    return finish();
}
