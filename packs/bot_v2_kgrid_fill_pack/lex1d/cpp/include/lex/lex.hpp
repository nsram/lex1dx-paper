// LEX: fast least-squares partition of 1-D data into k contiguous bins (scalar quantization).
// The Scala reference in ../scala implements the same algorithm; see docs/algorithm.md.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lex {

/// Which partition refinement started from.
enum class Start { Collapse, EqualFrequency, EqualWidth, Exact, CoarseExact };

/// Options for lex_plain and lex_hybrid (ckmeans_exact takes none).
struct Options {
    /// 0 groups exact duplicates only. A positive value joins each sorted value to the current
    /// run if it lies within `precision` of the run's mean (the 2024 deck rule).
    double precision = 0.0;
    /// Also try equal-frequency and equal-width starts; refine whichever has the lowest SSE.
    bool multi_start = false;
    /// A move is applied only if it lowers SSE by more than tolerance x (size of its terms).
    double tolerance = 1e-11;
    /// lex_hybrid only (lex_plain ignores it): Collapse stops at M = max(k, min(coarse_groups, atoms))
    /// groups, the exact DP cuts those groups into k bins, and refinement continues from there.
    /// M == k (coarse_groups <= k, including 0) is lex_plain's Collapse start; M == atoms starts
    /// refinement from the global optimum. With multi_start, this candidate replaces the Collapse
    /// candidate; equal-frequency and equal-width are still tried. The default is 4096: on
    /// continuous data up to N = 1e6 the relative excess SSE over exact stayed below 0.002% for
    /// k <= 64 and 0.08% for k <= 512; at N = 1e6 each doubling of M lowered it about 3.9x (README).
    std::size_t coarse_groups = 4096;
};

/// The 2024 deck setting (1e-4 grouping, three-way start). The deck preset is lex_plain(xs, k, deck()).
inline Options deck() { Options o; o.precision = 1e-4; o.multi_start = true; return o; }

/// A fitted partition of the sorted data. Bin i holds sorted values bounds[i] .. bounds[i+1]-1.
struct Result {
    std::vector<int> bounds;     ///< k + 1 cut indices into the sorted data; bounds[0] = 0, bounds[k] = n
    std::vector<int> sizes;      ///< values per bin
    std::vector<double> centers; ///< bin means (of the raw values)
    std::vector<double> upper;   ///< largest value in each bin
    double sse = 0;              ///< sum of squared deviations from the bin means (raw values)
    Start start = Start::Collapse;
    int atoms = 0;  ///< groups after duplicate / precision grouping
    int sweeps = 0; ///< refinement sweeps
    int moves = 0;  ///< refinement moves applied
    bool presorted = false; ///< the input was already in ascending order, so the sort was skipped

    int k() const { return static_cast<int>(sizes.size()); } ///< bins used: min(requested k, atoms)
    /// Bin index for x: the first bin whose upper value is >= x (the last bin if none).
    int assign(double x) const;
};

/// lex_hybrid: the coarse exact start (Collapse to options.coarse_groups groups, the exact DP on
/// those groups) followed by LEX refinement. This is the default method. With at most
/// coarse_groups atoms it returns ckmeans_exact's partition.
/// lex_plain: plain LEX, the Collapse start followed by refinement (options.coarse_groups ignored).
/// Both return a local optimum: no 3-bin-window move improves it. O(n log n) sort plus about
/// O(d log d + k m log m), m = coarse groups (lex_plain: m = k).
/// Already-sorted input is detected by an O(n) pass and the sort is skipped (Result::presorted).
///
/// Ties: among partitions with equal SSE, the one with the largest cuts (compared from the last
/// cut down) is returned, as Ckmeans.1d.dp (R) does. The DP scans split points right to left with
/// a strict <; after convergence Refine runs tie sweeps that shift bounds right when
/// |delta| <= 4 ulp x (size of its terms), exactly 0 with weights; ckmeans_exact runs the same pass.
/// Equal values always share an atom (-0.0 is +0.0); Collapse merges the leftmost pair on exactly
/// equal cost; with weights, equal values are ordered by weight. Results never depend on input order.
/// k larger than the number of atoms (distinct values, or precision groups) is clamped to it, as
/// Ckmeans.1d.dp does; Result::k() is the number of bins actually used.
/// Throws std::invalid_argument on empty input, non-finite values, k < 1 or bad options, and on
/// finite input too large for double arithmetic: n * max |x - mean| (with weights, max(W, 1) *
/// max |x - mean| and W) above sqrt(DBL_MAX) / 4 = 3.35e153. The message gives the limit and the
/// observed values (docs/algorithm.md, "Overflow").
Result lex_hybrid(const double* xs, std::size_t n, int k, const Options& options = {});
inline Result lex_hybrid(const std::vector<double>& xs, int k, const Options& options = {}) {
    return lex_hybrid(xs.data(), xs.size(), k, options);
}
Result lex_plain(const double* xs, std::size_t n, int k, const Options& options = {});
inline Result lex_plain(const std::vector<double>& xs, int k, const Options& options = {}) {
    return lex_plain(xs.data(), xs.size(), k, options);
}

/// Global optimum by dynamic programming (Fisher 1958) with divide and conquer, as Ckmeans.1d.dp:
/// O(k d log d) time and O(k d) memory for d distinct values.
Result ckmeans_exact(const double* xs, std::size_t n, int k);
inline Result ckmeans_exact(const std::vector<double>& xs, int k) { return ckmeans_exact(xs.data(), xs.size(), k); }

/// Weighted variants. `weights` is a non-negative finite weight per value (same length as `xs`).
/// Passing a null weight pointer, or using the overloads above, keeps the unweighted solver and
/// is bit-identical to previous releases. `Result::sse` is then the weighted SSE, and
/// `Result::centers` are the weighted bin means. A bin whose values all have weight 0 has no
/// weighted mean; its centre is the weighted mean of all the data, and it adds 0 to the SSE.
Result lex_hybrid(const double* xs, std::size_t n, int k, const double* weights, const Options& options);
inline Result lex_hybrid(const std::vector<double>& xs, const std::vector<double>& weights, int k, const Options& options = {}) {
    if (weights.size() != xs.size()) throw std::invalid_argument("weights size must match input");
    return lex_hybrid(xs.data(), xs.size(), k, weights.data(), options);
}
Result lex_plain(const double* xs, std::size_t n, int k, const double* weights, const Options& options);
inline Result lex_plain(const std::vector<double>& xs, const std::vector<double>& weights, int k, const Options& options = {}) {
    if (weights.size() != xs.size()) throw std::invalid_argument("weights size must match input");
    return lex_plain(xs.data(), xs.size(), k, weights.data(), options);
}
Result ckmeans_exact(const double* xs, std::size_t n, int k, const double* weights);
inline Result ckmeans_exact(const std::vector<double>& xs, const std::vector<double>& weights, int k) {
    if (weights.size() != xs.size()) throw std::invalid_argument("weights size must match input");
    return ckmeans_exact(xs.data(), xs.size(), k, weights.data());
}

} // namespace lex
