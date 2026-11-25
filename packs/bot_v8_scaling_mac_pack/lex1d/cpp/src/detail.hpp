// Internal pieces of LEX, one per source file. Not part of the public API.
#pragma once

#include <functional>
#include <vector>

namespace lex::detail {

/// Sorted data, centred on its mean and grouped into atoms: runs of values treated as one value.
/// Atom j covers elements start[j] .. start[j+1]-1, all with centred value value[j]; atom_id[e] is
/// the atom of element e. Bin sums come from double-double prefix sums, so a bin's mean depends
/// only on its bounds (not on the history of moves) and is accurate to a few ulps.
struct Atoms {
    std::vector<double> sorted;
    double center = 0;
    std::vector<int> atom_id, start;
    std::vector<double> value;
    std::vector<double> prefix_hi, prefix_lo;
    int n = 0, d = 0;
    bool presorted = false; ///< set by the caller: the input was already sorted (no sort ran)
    /// Weighted path only. The unweighted constructor leaves these at their defaults and is
    /// unchanged; mean / merge_cost / sse then use the original formulas. The weighted constructor
    /// leaves prefix_hi / prefix_lo empty (mean() reads w_* / ws_* instead), so code that reads
    /// prefix_hi directly must check `weighted` first.
    bool weighted = false;
    std::vector<double> elem_weight; // per sorted element; empty if unweighted
    std::vector<double> w_hi, w_lo, ws_hi, ws_lo, wq_hi, wq_lo;

    /// Groups sorted values: precision 0 joins exact duplicates; precision > 0 joins a value to
    /// the current run if it is within `precision` of the run's mean.
    Atoms(std::vector<double> sorted_values, double precision);
    /// Same with a non-negative weight per sorted value (aligned with `sorted_values`).
    Atoms(std::vector<double> sorted_values, std::vector<double> sorted_weights, double precision);

    /// Sum of weights / weighted centred sums of elements l .. r-1 (weighted path).
    double wsum(int l, int r) const {
        return (w_hi[r] - w_hi[l]) + (w_lo[r] - w_lo[l]);
    }
    double wS(int l, int r) const {
        return (ws_hi[r] - ws_hi[l]) + (ws_lo[r] - ws_lo[l]);
    }
    double wQ(int l, int r) const {
        return (wq_hi[r] - wq_hi[l]) + (wq_lo[r] - wq_lo[l]);
    }
    /// Weighted SSE of elements l .. r-1 about their weighted mean (0 if the weight is 0).
    double range_sse(int l, int r) const {
        const double W = wsum(l, r);
        if (W == 0) return 0;
        const double S = wS(l, r), Q = wQ(l, r);
        const double s = Q - S * S / W;
        return s > 0 ? s : 0;
    }

    /// Mean centred value of elements l .. r-1 (l < r). Weighted: if the range has total weight 0
    /// the mean is 0, i.e. the global weighted mean `center` (a bin of zero-weight values only).
    double mean(int l, int r) const {
        if (weighted) {
            const double W = wsum(l, r);
            return W == 0 ? 0.0 : wS(l, r) / W;
        }
        return ((prefix_hi[r] - prefix_hi[l]) + (prefix_lo[r] - prefix_lo[l])) / (r - l);
    }
    /// SSE increase from merging bins l..m-1 and m..r-1: na nb / (na + nb) (mean_a - mean_b)^2.
    double merge_cost(int l, int m, int r) const {
        if (weighted) {
            const double W1 = wsum(l, m), W2 = wsum(m, r);
            if (W1 == 0 || W2 == 0) return 0;
            const double diff = mean(l, m) - mean(m, r);
            return W1 * W2 / (W1 + W2) * diff * diff;
        }
        double na = m - l, nb = r - m;
        double diff = mean(l, m) - mean(m, r);
        return na * nb / (na + nb) * diff * diff;
    }
    /// SSE of the atom values for bounds b (two passes per bin).
    double sse(const std::vector<int>& b) const;
};

/// Greedy adjacent merge down to k bins; exactly equal costs merge the leftmost pair (this fixes
/// only the start; Refine's tie sweeps decide final ties). Returns bounds.
std::vector<int> collapse(const Atoms& a, int k);

/// 3-bin-window refinement of b in place. on_move (tests only) sees b after each applied move.
struct RefineStats { int sweeps = 0, moves = 0; };
/// Most gain/tie alternations in refine (ties drift to the largest cuts).
constexpr int kTieRounds = 8;
/// A tie move's computed |delta| is at most kTieEps * mag: equal SSE up to the rounding of the
/// unweighted delta (4 ulp of its terms). Weighted refinement counts only exactly zero deltas.
constexpr double kTieEps = 4 * 2.220446049250313e-16;
RefineStats refine(const Atoms& a, std::vector<int>& b, double tolerance,
                   const std::function<void(const std::vector<int>&)>& on_move = {});

/// Globally optimal bounds by DP with divide and conquer.
std::vector<int> exact(const Atoms& a, int k);
/// Same over weighted items: item i covers elements start[i] .. start[i+1]-1 and has centred mean
/// value[i]. Optimal among partitions that never split an item. Returns element bounds.
std::vector<int> exact(const std::vector<int>& start, const std::vector<double>& value, int k);
/// Weighted exact DP over the groups with element bounds `start`, using the Atoms prefixes.
/// Unweighted Atoms still go through exact(start, value, k).
std::vector<int> exact_weighted(const Atoms& a, const std::vector<int>& start, int k);

/// Canonicalizes -0.0 to +0.0, then sorts ascending unless the values are already sorted (an O(n)
/// check). Returns true if a sort was needed.
bool sort_values(std::vector<double>& v);
/// Same for (value, weight) pairs, ordered by value and then by weight. x and w have equal length.
bool sort_pairs(std::vector<double>& x, std::vector<double>& w);

std::vector<int> equal_frequency(int n, int k);
/// Bins of equal value width; empty if any bin would be empty. The last bin always ends at n.
std::vector<int> equal_width(const std::vector<double>& sorted, int k);

} // namespace lex::detail
