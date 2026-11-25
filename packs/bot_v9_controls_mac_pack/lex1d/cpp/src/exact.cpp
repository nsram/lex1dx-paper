// Exact 1-D least-squares partition (optimal k-means in one dimension) by dynamic programming over
// atoms (runs of equal values), as in Fisher (1958). Each layer is filled by divide and conquer on
// the monotone optimal split point, which holds because the bin SSE cost satisfies the quadrangle
// inequality (Wu 1991; Grønlund et al. 2017).
//
// Time O(k d log d) after the sort; memory O(k d) ints for back-pointers plus O(d) doubles.
// Bin costs use double-double prefix sums of centred values and squares: Q - S^2/n.
// Ties: among exactly equal costs the largest split index wins (as Ckmeans.1d.dp does), so the
// backtrack returns the optimum whose cuts are largest, compared from the last cut down.
#include "detail.hpp"

#include <algorithm>
#include <limits>

namespace lex::detail {
namespace {

void two_sum(std::vector<double>& hi, std::vector<double>& lo, int j, double x) {
    const double s = hi[j] + x;
    const double bv = s - hi[j];
    hi[j + 1] = s;
    lo[j + 1] = lo[j] + (hi[j] - (s - bv)) + (x - bv);
}

struct Dp {
    Dp(const std::vector<int>& st, const std::vector<double>& v) : start(st), value(v) {}
    const std::vector<int>& start;   // item i covers elements start[i] .. start[i+1]-1
    const std::vector<double>& value; // centred value of item i
    std::vector<double> s_hi, s_lo, q_hi, q_lo; // prefix sums over atoms of count*value, count*value^2
    std::vector<double> prev, cur;               // best SSE of atoms 0..j-1 in b-1 and b bins
    std::vector<int>* from = nullptr;            // row of back-pointers for the current layer

    /// SSE of items i .. j-1 (i < j).
    double cost(int i, int j) const {
        const double n = start[j] - start[i];
        const double s = (s_hi[j] - s_hi[i]) + (s_lo[j] - s_lo[i]);
        const double q = (q_hi[j] - q_hi[i]) + (q_lo[j] - q_lo[i]);
        return std::max(0.0, q - s * s / n);
    }

    /// Fill cur[j] for j in [lo, hi], knowing the optimal split lies in [opt_lo, opt_hi].
    void fill(int lo, int hi, int opt_lo, int opt_hi) {
        if (lo > hi) return;
        const int mid = (lo + hi) >> 1;
        double best = std::numeric_limits<double>::infinity();
        int arg = std::min(opt_hi, mid - 1);
        const int last = std::min(opt_hi, mid - 1);
        // right to left, strict <: among exactly equal costs the largest split index wins
        for (int i = last; i >= opt_lo; --i) {
            const double c = prev[i] + cost(i, mid);
            if (c < best) { best = c; arg = i; }
        }
        cur[mid] = best;
        (*from)[mid] = arg;
        fill(lo, mid - 1, opt_lo, arg);
        fill(mid + 1, hi, arg, opt_hi);
    }
};

} // namespace

std::vector<int> exact(const std::vector<int>& start, const std::vector<double>& value, int k) {
    const int d = static_cast<int>(value.size());
    Dp dp(start, value);
    for (auto* v : {&dp.s_hi, &dp.s_lo, &dp.q_hi, &dp.q_lo}) v->assign(d + 1, 0.0);
    for (int j = 0; j < d; ++j) {
        const double c = start[j + 1] - start[j];
        const double v = value[j];
        two_sum(dp.s_hi, dp.s_lo, j, c * v);
        two_sum(dp.q_hi, dp.q_lo, j, c * v * v);
    }
    dp.prev.assign(d + 1, 0.0);
    for (int j = 1; j <= d; ++j) dp.prev[j] = dp.cost(0, j); // 1 bin
    dp.cur.assign(d + 1, 0.0);

    // from[b-2][j]: start item of the last of b bins covering atoms 0..j-1 (b = 2..k)
    std::vector<std::vector<int>> from(k > 1 ? k - 1 : 0, std::vector<int>(d + 1));
    for (int b = 2; b <= k; ++b) {
        dp.from = &from[b - 2];
        dp.fill(b, d, b - 1, d - 1);
        std::swap(dp.prev, dp.cur);
    }

    std::vector<int> bounds(k + 1, 0);
    bounds[k] = start[d];
    int end = d;
    for (int b = k; b >= 2; --b) {
        end = from[b - 2][end];
        bounds[b - 1] = start[end];
    }
    return bounds;
}

std::vector<int> exact(const Atoms& a, int k) {
    if (a.weighted) return exact_weighted(a, a.start, k);
    return exact(a.start, a.value, k);
}

std::vector<int> exact_weighted(const Atoms& a, const std::vector<int>& start, int k) {
    const int d = static_cast<int>(start.size()) - 1;
    // the prefix values at the d+1 group bounds, copied once into contiguous arrays (the same
    // doubles, so the same costs as a.range_sse(start[i], start[j]))
    std::vector<double> wh(d + 1), wl(d + 1), sh(d + 1), sl(d + 1), qh(d + 1), ql(d + 1);
    for (int j = 0; j <= d; ++j) {
        const int e = start[j];
        wh[j] = a.w_hi[e]; wl[j] = a.w_lo[e]; sh[j] = a.ws_hi[e]; sl[j] = a.ws_lo[e]; qh[j] = a.wq_hi[e]; ql[j] = a.wq_lo[e];
    }
    auto cost = [&](int i, int j) { // = a.range_sse(start[i], start[j])
        const double W = (wh[j] - wh[i]) + (wl[j] - wl[i]);
        if (W == 0) return 0.0;
        const double S = (sh[j] - sh[i]) + (sl[j] - sl[i]), Q = (qh[j] - qh[i]) + (ql[j] - ql[i]);
        const double s = Q - S * S / W;
        return s > 0 ? s : 0.0;
    };
    std::vector<double> prev(d + 1, 0.0), cur(d + 1, 0.0);
    for (int j = 1; j <= d; ++j) prev[j] = cost(0, j);
    const std::size_t row = static_cast<std::size_t>(d) + 1;
    std::vector<int> from(k > 1 ? (k - 1) * row : 0);
    auto fill = [&](auto& self, int lo, int hi, int opt_lo, int opt_hi, int* fr) -> void {
        if (lo > hi) return;
        const int mid = (lo + hi) >> 1;
        double best = std::numeric_limits<double>::infinity();
        const int last = std::min(opt_hi, mid - 1);
        int arg = last;
        for (int i = last; i >= opt_lo; --i) {
            const double c = prev[i] + cost(i, mid);
            if (c < best) { best = c; arg = i; }
        }
        cur[mid] = best;
        fr[mid] = arg;
        self(self, lo, mid - 1, opt_lo, arg, fr);
        self(self, mid + 1, hi, arg, opt_hi, fr);
    };
    for (int b = 2; b <= k; ++b) {
        fill(fill, b, d, b - 1, d - 1, from.data() + (b - 2) * row);
        std::swap(prev, cur);
    }
    std::vector<int> bounds(k + 1, 0);
    bounds[k] = start[d];
    int end = d;
    for (int b = k; b >= 2; --b) {
        end = from[(b - 2) * row + end];
        bounds[b - 1] = start[end];
    }
    return bounds;
}

} // namespace lex::detail
