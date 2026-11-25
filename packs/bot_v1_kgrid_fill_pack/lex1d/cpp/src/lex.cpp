#include "lex/lex.hpp"
#include "detail.hpp"

#include <algorithm>
#include <charconv>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace lex {
namespace {

using detail::Atoms;

/// Java's Double.toString spelling for the values that appear in error messages.
std::string show(double x) {
    if (std::isnan(x)) return "NaN";
    if (std::isinf(x)) return x > 0 ? "Infinity" : "-Infinity";
    char buf[32];
    auto res = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, res.ptr);
}

[[noreturn]] void fail(const std::string& msg) { throw std::invalid_argument(msg); }

/// The overflow rule (docs/algorithm.md, "Overflow"), checked once on the input before any solver
/// runs. Let Z = max |x - mean| (the mean the solvers centre on) and N = n, or the total weight W.
/// If max(N, 1) * Z <= sqrt(DBL_MAX) / 4 and N <= sqrt(DBL_MAX) / 4, the sums the solvers form stay
/// far below DBL_MAX. A mean that is not finite (the sum of the values overflows) gives Z = inf and
/// fails too. Input that passes runs exactly as before; nothing inside the solvers checks again.
void check_range(const Atoms& a) {
    const double limit = std::sqrt(std::numeric_limits<double>::max()) / 4; // 3.351951982485649e+153
    const double z = std::max(std::abs(a.sorted.front() - a.center), std::abs(a.sorted.back() - a.center));
    const double total = a.weighted ? a.wsum(0, a.n) : static_cast<double>(a.n);
    const double product = std::max(total, 1.0) * z;
    if (product <= limit && total <= limit) return;
    const std::string lim = "the limit sqrt(DBL_MAX) / 4 = " + show(limit);
    const std::string zs = "max |x - mean| = " + show(z);
    std::string msg = "input range too large for double arithmetic: ";
    if (!(total <= limit))
        msg += "the total weight W = " + show(total) + " exceeds " + lim + " (" + zs + ")";
    else if (a.weighted)
        msg += "max(W, 1) * max |x - mean| = " + show(product) + " exceeds " + lim + " (W = " + show(total) + ", " + zs + ")";
    else
        msg += "n * max |x - mean| = " + show(product) + " exceeds " + lim + " (n = " + std::to_string(a.n) + ", " + zs + ")";
    if (!std::isfinite(a.center))
        msg += "; the mean is not finite because the sum of the values overflows (limit " +
               show(std::numeric_limits<double>::max()) + ")";
    fail(msg);
}

/// Validates input at the boundary and returns sorted, grouped data. Clamps k to the number of atoms.
Atoms prepare(const double* xs, std::size_t n, int& k, const Options& o) {
    if (n == 0) fail("input is empty");
    if (n >= static_cast<std::size_t>(INT_MAX)) fail("input has " + std::to_string(n) + " values; at most 2147483646 are supported");
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(xs[i]))
            fail("input value at index " + std::to_string(i) + " is " + show(xs[i]) + "; only finite values are allowed");
    if (k < 1) fail("k must be at least 1, got " + std::to_string(k));
    if (!(o.precision >= 0 && std::isfinite(o.precision))) fail("precision must be finite and >= 0, got " + show(o.precision));
    if (!(o.tolerance > 0 && o.tolerance < 1)) fail("tolerance must be in (0, 1), got " + show(o.tolerance));

    std::vector<double> sorted(xs, xs + n);
    const bool did_sort = detail::sort_values(sorted); // skipped if already sorted
    Atoms atoms(std::move(sorted), o.precision);
    atoms.presorted = !did_sort;
    check_range(atoms);
    if (k > atoms.d) k = atoms.d; // more bins than atoms: clamp, as Ckmeans.1d.dp does (Result::k() reports it)
    return atoms;
}

Result summarize(const Atoms& a, std::vector<int> b, Start start, detail::RefineStats stats) {
    const int k = static_cast<int>(b.size()) - 1;
    Result r;
    r.sizes.resize(k);
    r.centers.resize(k);
    r.upper.resize(k);
    double sse = 0;
    for (int i = 0; i < k; ++i) { // per bin: two passes over its values, centred
        const int lo = b[i], hi = b[i + 1];
        double s = 0;
        for (int e = lo; e < hi; ++e) s += a.sorted[e] - a.center;
        const double m = s / (hi - lo);
        for (int e = lo; e < hi; ++e) {
            const double dv = a.sorted[e] - a.center - m;
            sse += dv * dv;
        }
        r.sizes[i] = hi - lo;
        r.centers[i] = a.center + m;
        r.upper[i] = a.sorted[hi - 1];
    }
    r.bounds = std::move(b);
    r.sse = sse;
    r.start = start;
    r.atoms = a.d;
    r.sweeps = stats.sweeps;
    r.moves = stats.moves;
    r.presorted = a.presorted;
    return r;
}

} // namespace

namespace detail {

std::vector<int> equal_frequency(int n, int k) { // first k - n%k bins get n/k values, the rest one more
    const int f = n / k, extra = n - f * k;
    std::vector<int> b(k + 1, 0);
    for (int i = 1; i <= k; ++i) b[i] = b[i - 1] + f + (i > k - extra ? 1 : 0);
    return b;
}

std::vector<int> equal_width(const std::vector<double>& sorted, int k) {
    const int n = static_cast<int>(sorted.size());
    const double lo = sorted[0];
    const double width = (sorted[n - 1] - lo) / k;
    std::vector<int> b(k + 1, 0);
    for (int i = 1; i < k; ++i) // number of values <= the cut
        b[i] = static_cast<int>(std::upper_bound(sorted.begin(), sorted.end(), lo + i * width) - sorted.begin());
    b[k] = n;
    for (int i = 1; i <= k; ++i)
        if (b[i] <= b[i - 1]) return {};
    return b;
}

} // namespace detail

Result lex_hybrid(const double* xs, std::size_t n, int k, const Options& options) {
    const Atoms a = prepare(xs, n, k, options);
    const std::size_t m = std::max<std::size_t>(k, std::min<std::size_t>(options.coarse_groups, a.d));
    std::vector<int> b;
    Start start = Start::Collapse;
    if (m <= static_cast<std::size_t>(k)) {
        b = detail::collapse(a, k);
    } else if (m == static_cast<std::size_t>(a.d)) { // every atom is a group: DP over the atoms
        b = detail::exact(a, k);
        start = Start::CoarseExact;
    } else { // Collapse to m groups, then the DP over the groups (O(m) means from the prefix sums)
        const std::vector<int> groups = detail::collapse(a, static_cast<int>(m));
        std::vector<double> means(m);
        for (std::size_t g = 0; g < m; ++g) means[g] = a.mean(groups[g], groups[g + 1]);
        b = detail::exact(groups, means, k);
        start = Start::CoarseExact;
    }
    if (options.multi_start) { // lowest atom-value SSE wins; ties keep the earlier candidate
        double best = a.sse(b);
        std::vector<int> ef = detail::equal_frequency(a.n, k);
        if (const double s = a.sse(ef); s < best) { best = s; b = std::move(ef); start = Start::EqualFrequency; }
        std::vector<int> ew = detail::equal_width(a.sorted, k);
        if (!ew.empty())
            if (const double s = a.sse(ew); s < best) { b = std::move(ew); start = Start::EqualWidth; }
    }
    const auto stats = detail::refine(a, b, options.tolerance);
    return summarize(a, std::move(b), start, stats);
}

namespace {
/// The DP's optimum, then the refine pass with its tie sweeps, so exact ties that the DP's rounding
/// hides still end on the largest cuts (as Ckmeans.1d.dp). Gain moves cannot apply at the optimum;
/// Result::sweeps stays 0 and Result::moves counts tie moves.
std::vector<int> exact_largest_cuts(const Atoms& a, int k, detail::RefineStats& stats) {
    std::vector<int> b = detail::exact(a, k);
    stats = detail::refine(a, b, Options{}.tolerance);
    stats.sweeps = 0;
    return b;
}
} // namespace

Result ckmeans_exact(const double* xs, std::size_t n, int k) {
    const Atoms a = prepare(xs, n, k, Options{});
    detail::RefineStats stats;
    std::vector<int> b = exact_largest_cuts(a, k, stats);
    return summarize(a, std::move(b), Start::Exact, stats);
}

int Result::assign(double x) const {
    const auto it = std::lower_bound(upper.begin(), upper.end(), x); // first upper >= x
    return it == upper.end() ? k() - 1 : static_cast<int>(it - upper.begin());
}

namespace {

Atoms prepare_weighted(const double* xs, const double* ws, std::size_t n, int& k, const Options& o) {
    if (n == 0) fail("input is empty");
    if (n >= static_cast<std::size_t>(INT_MAX)) fail("input has " + std::to_string(n) + " values; at most 2147483646 are supported");
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(xs[i]))
            fail("input value at index " + std::to_string(i) + " is " + show(xs[i]) + "; only finite values are allowed");
    if (ws == nullptr) fail("weights pointer is null");
    double W = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(ws[i]) || ws[i] < 0)
            fail("weight at index " + std::to_string(i) + " is " + show(ws[i]) + "; weights must be finite and >= 0");
        W += ws[i];
    }
    if (!(W > 0) || !std::isfinite(W)) fail("weights must sum to a positive finite total, got " + show(W));
    if (k < 1) fail("k must be at least 1, got " + std::to_string(k));
    if (!(o.precision >= 0 && std::isfinite(o.precision))) fail("precision must be finite and >= 0, got " + show(o.precision));
    if (!(o.tolerance > 0 && o.tolerance < 1)) fail("tolerance must be in (0, 1), got " + show(o.tolerance));

    std::vector<double> sorted(xs, xs + n), weights(ws, ws + n);
    const bool did_sort = detail::sort_pairs(sorted, weights); // skipped if already sorted
    Atoms atoms(std::move(sorted), std::move(weights), o.precision);
    atoms.presorted = !did_sort;
    check_range(atoms);
    if (k > atoms.d) k = atoms.d; // clamp, as in prepare()
    return atoms;
}

Result summarize_weighted(const Atoms& a, std::vector<int> b, Start start, detail::RefineStats stats) {
    const int k = static_cast<int>(b.size()) - 1;
    Result r;
    r.sizes.resize(k);
    r.centers.resize(k);
    r.upper.resize(k);
    double sse = 0;
    for (int i = 0; i < k; ++i) {
        const int lo = b[i], hi = b[i + 1];
        double W = 0, S = 0;
        for (int e = lo; e < hi; ++e) { W += a.elem_weight[e]; S += a.elem_weight[e] * (a.sorted[e] - a.center); }
        // a bin whose values all have weight 0 has no weighted mean: its centre falls back to the
        // global weighted mean (centred 0) and it adds nothing to the SSE
        const double m = W == 0 ? 0.0 : S / W;
        if (W > 0) {
            for (int e = lo; e < hi; ++e) {
                const double dv = a.sorted[e] - a.center - m;
                sse += a.elem_weight[e] * dv * dv;
            }
        }
        r.sizes[i] = hi - lo;
        r.centers[i] = a.center + m;
        r.upper[i] = a.sorted[hi - 1];
    }
    r.bounds = std::move(b);
    r.sse = sse;
    r.start = start;
    r.atoms = a.d;
    r.sweeps = stats.sweeps;
    r.moves = stats.moves;
    r.presorted = a.presorted;
    return r;
}

} // namespace

Result lex_hybrid(const double* xs, std::size_t n, int k, const double* weights, const Options& options) {
    if (weights == nullptr) return lex_hybrid(xs, n, k, options); // unweighted path, unchanged
    const Atoms a = prepare_weighted(xs, weights, n, k, options);
    const std::size_t m = std::max<std::size_t>(k, std::min<std::size_t>(options.coarse_groups, a.d));
    std::vector<int> b;
    Start start = Start::Collapse;
    if (m <= static_cast<std::size_t>(k)) {
        b = detail::collapse(a, k);
    } else if (m == static_cast<std::size_t>(a.d)) {
        b = detail::exact(a, k);
        start = Start::CoarseExact;
    } else {
        const std::vector<int> groups = detail::collapse(a, static_cast<int>(m));
        b = detail::exact_weighted(a, groups, k);
        start = Start::CoarseExact;
    }
    if (options.multi_start) {
        double best = a.sse(b);
        std::vector<int> ef = detail::equal_frequency(a.n, k);
        if (const double s = a.sse(ef); s < best) { best = s; b = std::move(ef); start = Start::EqualFrequency; }
        std::vector<int> ew = detail::equal_width(a.sorted, k);
        if (!ew.empty())
            if (const double s = a.sse(ew); s < best) { b = std::move(ew); start = Start::EqualWidth; }
    }
    const auto stats = detail::refine(a, b, options.tolerance);
    return summarize_weighted(a, std::move(b), start, stats);
}

Result ckmeans_exact(const double* xs, std::size_t n, int k, const double* weights) {
    if (weights == nullptr) return ckmeans_exact(xs, n, k);
    const Atoms a = prepare_weighted(xs, weights, n, k, Options{});
    detail::RefineStats stats;
    std::vector<int> b = exact_largest_cuts(a, k, stats);
    return summarize_weighted(a, std::move(b), Start::Exact, stats);
}

namespace {
Options without_coarse(Options o) { o.coarse_groups = 0; return o; }
} // namespace

Result lex_plain(const double* xs, std::size_t n, int k, const Options& options) {
    return lex_hybrid(xs, n, k, without_coarse(options));
}

Result lex_plain(const double* xs, std::size_t n, int k, const double* weights, const Options& options) {
    return lex_hybrid(xs, n, k, weights, without_coarse(options));
}

} // namespace lex
