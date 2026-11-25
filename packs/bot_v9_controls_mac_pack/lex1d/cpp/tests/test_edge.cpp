#include "support.hpp"
#include "../src/detail.hpp"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

using namespace lex;

static void rejects(const std::vector<double>& xs, int k, const std::string& fragment, Options o = {}) {
    try {
        solve(xs, k, o);
        CHECK_MSG(false, "no exception for " + fragment);
    } catch (const std::invalid_argument& e) {
        CHECK_MSG(std::string(e.what()).find(fragment) != std::string::npos, e.what());
    }
}

int main() {
    const double inf = std::numeric_limits<double>::infinity();
    test("invalid input is rejected with a clear message", [&] {
        rejects({}, 1, "empty");
        rejects({1.0, std::nan("")}, 1, "index 1 is NaN");
        rejects({inf}, 1, "index 0 is Infinity");
        rejects({0.0, -inf}, 1, "index 1 is -Infinity");
        rejects({1.0, 2.0}, 0, "k must be at least 1");
        rejects({1.0}, 1, "precision", Options{-1});
        rejects({1.0}, 1, "precision", Options{inf});
        rejects({1.0}, 1, "tolerance", Options{0, false, 0});
        try { ckmeans_exact(std::vector<double>{std::nan("")}, 1); CHECK(false); } catch (const std::invalid_argument&) {}
    });

    test("input too large for double arithmetic is rejected: n * max |x - mean| <= sqrt(DBL_MAX) / 4", [&] {
        const double L = std::sqrt(DBL_MAX) / 4;
        CHECK(L == 3.351951982485649e153);
        const std::string lim = "exceeds the limit sqrt(DBL_MAX) / 4 = 3.351951982485649e+153";
        auto all = [&](const std::vector<double>& xs, int k, const std::string& fragment) {
            rejects(xs, k, fragment);
            rejects(xs, k, fragment, plain());
            rejects(xs, k, fragment, deck());
            try { ckmeans_exact(xs, k); CHECK_MSG(false, "no exception for " + fragment); }
            catch (const std::invalid_argument& e) { CHECK_MSG(std::string(e.what()).find(fragment) != std::string::npos, e.what()); }
        };
        auto passes = [&](const std::vector<double>& xs, int k) {
            for (const Options& o : {Options{}, plain(), deck()}) CHECK(std::isfinite(solve(xs, k, o).sse));
            CHECK(std::isfinite(ckmeans_exact(xs, k).sse));
        };
        // n = 2 at the limit: {-a, a} has mean exactly 0, so n * max |x - mean| = 2a
        const double a = L / 2, a_up = std::nextafter(a, DBL_MAX);
        passes({-a, a}, 2);
        all({-a_up, a_up}, 2, "n * max |x - mean| = 3.3519519824856493e+153 " + lim);
        all({-a_up, a_up}, 2, "(n = 2, max |x - mean| = 1.6759759912428246e+153)");
        // 1000 values: half at -2^505, half at +2^505 (about 1.05e152; the sums are exact, so the mean is 0 on any platform)
        std::vector<double> runs;
        for (int i = 0; i < 1000; ++i) runs.push_back(i < 500 ? -std::ldexp(1.0, 505) : std::ldexp(1.0, 505));
        all(runs, 2, lim + " (n = 1000, max |x - mean| = 1.0474849945267654e+152)");
        // a mean that is not finite: the sum of the values overflows
        all({1e308, 1e308}, 1, "the mean is not finite because the sum of the values overflows (limit 1.7976931348623157e+308)");
        all({1.7e308, -1e308, 1.7e308, -1e308}, 2, "mean is not finite");
        // the case that broke Refine's rounding bound before (cuts 3|1, SSE 1.04e308) is above the limit
        const double b = std::sqrt(DBL_MAX) * 0.93 / 2;
        all({-b, -b, b, b * (1 + 1.0 / 256)}, 2, lim);
        // data just below the limit match the same data scaled by 2^-600 (exact scaling)
        const double c = L / 1000 / 1.2; // mean c / 6, so n * max |x - mean| = 0.97 L
        std::vector<double> edge, scaled;
        for (int i = 0; i < 1000; ++i) edge.push_back(i % 3 == 0 ? -c : i % 3 == 1 ? c * 0.5 : c);
        for (double x : edge) scaled.push_back(std::ldexp(x, -600));
        for (int k : {2, 3, 5}) {
            for (const Options& o : {Options{}, plain()}) { // not deck: its precision does not scale
                const Result r = solve(edge, k, o), q = solve(scaled, k, o);
                CHECK(r.bounds == q.bounds && r.moves == q.moves && r.sweeps == q.sweeps && std::isfinite(r.sse));
            }
            CHECK(ckmeans_exact(edge, k).bounds == ckmeans_exact(scaled, k).bounds);
        }
        // weighted: max(W, 1) * max |x - mean| and W against the same limit
        auto wrejects = [&](const std::vector<double>& xs, const std::vector<double>& ws, const std::string& fragment) {
            for (int m = 0; m < 3; ++m) {
                try {
                    if (m == 0) lex_hybrid(xs, ws, 2); else if (m == 1) lex_plain(xs, ws, 2); else ckmeans_exact(xs, ws, 2);
                    CHECK_MSG(false, "no exception for " + fragment);
                } catch (const std::invalid_argument& e) { CHECK_MSG(std::string(e.what()).find(fragment) != std::string::npos, e.what()); }
            }
        };
        auto wpasses = [&](const std::vector<double>& xs, const std::vector<double>& ws) {
            CHECK(std::isfinite(lex_hybrid(xs, ws, 2).sse) && std::isfinite(lex_plain(xs, ws, 2).sse) &&
                  std::isfinite(ckmeans_exact(xs, ws, 2).sse));
        };
        const double a4 = L / 8; // W = 8: 8 * a4 = L
        wpasses({-a4, a4}, {4, 4});
        wrejects({-std::nextafter(a4, DBL_MAX), std::nextafter(a4, DBL_MAX)}, {4, 4}, "max(W, 1) * max |x - mean| = 3.3519519824856493e+153 " + lim + " (W = 8, ");
        wpasses({-L, L}, {0.25, 0.25}); // W < 1 counts as 1
        wrejects({-std::nextafter(L, DBL_MAX), std::nextafter(L, DBL_MAX)}, {0.25, 0.25}, "(W = 0.5, max |x - mean| = 3.3519519824856493e+153)");
        wrejects({0, 1}, {L, L}, "the total weight W = 6.703903964971298e+153 " + lim + " (max |x - mean| = 0.5)");
        wrejects({-1e200, 0, 1e200}, {1e10, 1, 1e10}, "max(W, 1) * max |x - mean| = 2.0000000001e+210 " + lim);
    });

    test("built without FMA contraction: a case whose SSE bits change with FMA", [] {
        // fixtures/data/fma_60.txt (random.Random(29).random() x 60). With -ffp-contract=fast
        // -march=haswell (x86-64, -O2 or -O3) the SSE comes out ...302 or the bounds move.
        const std::vector<double> xs{
        0x1.18a30f84df15ap-1, 0x1.6221bccc9caacp-2, 0x1.b0905266c46fcp-1, 0x1.2786160b1a082p-2,
        0x1.054bf17b583cdp-1, 0x1.6010d9b24b0e0p-2, 0x1.a9773d826331ap-2, 0x1.f29d9ec211ecdp-1,
        0x1.a8e39bb1e1138p-4, 0x1.c7604c537b312p-2, 0x1.cc67e8cd9645cp-3, 0x1.65f4c70af257ep-2,
        0x1.ff5add1fdd6dfp-1, 0x1.5133bf756f242p-2, 0x1.372ea662b485dp-1, 0x1.adb00b6b04abap-2,
        0x1.75fe044cc203cp-1, 0x1.c9f0c47d11392p-2, 0x1.5d0d70382db68p-3, 0x1.bc539be45e6dep-2,
        0x1.a1da83564e994p-2, 0x1.95df146bc3252p-1, 0x1.200df3e181022p-1, 0x1.923ad58b2b590p-3,
        0x1.59023d5b091a4p-1, 0x1.318953239e854p-2, 0x1.06168d7007fccp-1, 0x1.c243ffd1999ecp-3,
        0x1.a49cd37ab4430p-2, 0x1.927d6bde3d622p-1, 0x1.080a5509320c5p-1, 0x1.813076d3a74b1p-1,
        0x1.83735d9327d79p-1, 0x1.2c9783ad4f515p-1, 0x1.4a52749e2b71bp-1, 0x1.96d2566915452p-2,
        0x1.3c6434b64ce97p-1, 0x1.8df231794c206p-1, 0x1.96419f764dab8p-4, 0x1.e0110638f0106p-1,
        0x1.780d176dc138cp-2, 0x1.1ce6ab439aff4p-3, 0x1.333350b045f6cp-3, 0x1.5cda2866a66acp-2,
        0x1.1bb6911aadc08p-4, 0x1.9768a094e7712p-2, 0x1.a90cfb5005fcap-2, 0x1.4f89939b439edp-1,
        0x1.ddc9312400a4ap-2, 0x1.deddcba784ba0p-1, 0x1.51999589f6e8cp-1, 0x1.48b9975f6be50p-2,
        0x1.bff267360b543p-1, 0x1.65e9e37e4abc8p-1, 0x1.ef4447a8a4467p-1, 0x1.77fbdeaec5d72p-1,
        0x1.ee664deb28679p-1, 0x1.a44964abf04acp-3, 0x1.1d0a5181f293ep-2, 0x1.f1d23e7052a88p-1,
        };
        const std::vector<int> want{0, 11, 34, 51, 60};
        for (int m = 0; m < 4; ++m) {
            const Result r = m == 0 ? lex_hybrid(xs, 4) : m == 1 ? lex_plain(xs, 4) : m == 2 ? lex_plain(xs, 4, deck()) : ckmeans_exact(xs, 4);
            std::uint64_t bits;
            std::memcpy(&bits, &r.sse, 8);
            CHECK_MSG(bits == 0x3fd083e27c924301ULL && r.bounds == want, "method " + std::to_string(m));
        }
    });

    test("weighted Atoms leave the unweighted prefix sums empty", [] {
        const detail::Atoms w({1.0, 2.0, 2.0}, {1.0, 0.5, 2.0}, 0);
        CHECK(w.weighted && w.prefix_hi.empty() && w.prefix_lo.empty() && w.w_hi.size() == 4);
        const detail::Atoms u({1.0, 2.0, 2.0}, 0);
        CHECK(!u.weighted && u.prefix_hi.size() == 4);
    });

    test("k above the number of atoms is clamped to it (as Ckmeans.1d.dp); k() reports the bins used", [] {
        auto same_fit = [](const Result& a, const Result& b) {
            return a.bounds == b.bounds && a.sizes == b.sizes && a.sse == b.sse && a.centers == b.centers &&
                   a.start == b.start && a.atoms == b.atoms && a.sweeps == b.sweeps && a.moves == b.moves;
        };
        const std::vector<double> xs{3.0, 1.0, 1.0, 2.0, 3.0, 3.0, 2.0}; // 3 distinct values
        for (const Options& o : {Options{}, plain(), plain_deck()}) {
            for (int k : {4, 7, 1000}) {
                const Result r = solve(xs, k, o);
                CHECK(r.k() == 3); CHECK(r.sizes == (std::vector<int>{2, 2, 3})); CHECK(r.sse == 0.0);
                CHECK(same_fit(r, solve(xs, 3, o)));
            }
        }
        const Result e = ckmeans_exact(xs, 10);
        CHECK(e.k() == 3); CHECK(same_fit(e, ckmeans_exact(xs, 3)));
        const Result one = ckmeans_exact(std::vector<double>{1.0, 1.0}, 2);
        CHECK(one.k() == 1); CHECK(one.bounds == (std::vector<int>{0, 2}));
        // precision groups count as atoms: 1.0 and 1.05 share a group at precision 0.1
        const Result g = solve(std::vector<double>{1.0, 1.05, 2.0}, 3, Options{0.1});
        CHECK(g.k() == 2); CHECK(g.atoms == 2); CHECK(g.sizes == (std::vector<int>{2, 1}));
        // weighted
        const std::vector<double> w{1, 2, 1, 2, 1, 1, 1};
        const Result rw = solve(xs, w, 5), rw3 = solve(xs, w, 3);
        CHECK(rw.k() == 3); CHECK(same_fit(rw, rw3));
        const Result ew = ckmeans_exact(xs, w, 5);
        CHECK(ew.k() == 3); CHECK(same_fit(ew, ckmeans_exact(xs, w, 3)));
    });

    test("all-equal values: k = 1 has SSE 0", [] {
        const Result r = solve(std::vector<double>(7, 3.5), 1);
        CHECK(r.sizes == std::vector<int>{7}); CHECK(r.sse == 0.0); CHECK(r.centers[0] == 3.5);
        CHECK(r.bounds == (std::vector<int>{0, 7})); CHECK(r.atoms == 1);
        const Result e = ckmeans_exact(std::vector<double>(7, 3.5), 1);
        CHECK(e.sse == 0.0); CHECK(e.start == Start::Exact);
    });

    test("k = 1 gives the total sum of squares", [] {
        const std::vector<double> xs{1, 2, 4, 8};
        CHECK_NEAR(solve(xs, 1).sse, sse_of(xs, {4}), 1e-12);
    });

    test("k = number of distinct values gives SSE 0 and keeps duplicates together", [] {
        const Result r = solve(std::vector<double>{3, 1, 2, 2, 5, 1}, 4);
        CHECK(r.sizes == (std::vector<int>{2, 2, 1, 1})); CHECK(r.sse == 0.0);
        CHECK(r.bounds == (std::vector<int>{0, 2, 4, 5, 6}));
    });

    test("tiny inputs, k = 2 refinement, and signed zeros", [] {
        CHECK(solve(std::vector<double>{0, 1}, 2).sizes == (std::vector<int>{1, 1}));
        CHECK(solve(std::vector<double>{0, 0, 1, 2}, 2).sizes == (std::vector<int>{2, 2})); // optimum
        CHECK(solve(std::vector<double>{-0.0, 0.0, 0.0}, 1).atoms == 1);
        CHECK(solve(std::vector<double>{42.0}, 1).sse == 0.0);
    });

    test("input is not modified; assign() maps every value to its bin", [] {
        const std::vector<double> xs{9, 1, 5, 2, 7, 3};
        const std::vector<double> copy = xs;
        const Result r = solve(xs, 3);
        CHECK(xs == copy);
        std::vector<double> s = xs;
        std::sort(s.begin(), s.end());
        int e = 0;
        for (int i = 0; i < r.k(); ++i)
            for (int c = 0; c < r.sizes[i]; ++c) CHECK(r.assign(s[e++]) == i);
        CHECK(r.assign(100.0) == r.k() - 1);
    });

    test("ties merge the leftmost pair first", [] {
        std::vector<double> xs(12);
        for (int i = 0; i < 12; ++i) xs[i] = i; // every adjacent merge costs the same
        const detail::Atoms a(xs, 0);
        CHECK(detail::collapse(a, 11) == (std::vector<int>{0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}));
    });

    test("deck preset and output invariants", [] {
        std::vector<double> xs;
        for (int i = 0; i < 1000; ++i) xs.push_back(std::sin(i * 0.37) * 10);
        for (const Options& o : {Options{}, plain(), plain_deck()}) {
            const Result r = solve(xs, 16, o);
            CHECK(r.bounds.size() == 17u); CHECK(r.bounds.front() == 0); CHECK(r.bounds.back() == 1000);
            int total = 0;
            for (int i = 0; i < 16; ++i) { CHECK(r.sizes[i] == r.bounds[i + 1] - r.bounds[i]); CHECK(r.sizes[i] > 0); total += r.sizes[i]; }
            CHECK(total == 1000);
            CHECK(std::is_sorted(r.upper.begin(), r.upper.end()));
        }
    });

    return finish();
}
