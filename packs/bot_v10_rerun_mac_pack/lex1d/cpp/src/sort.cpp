// Sorting with an O(n) sorted-input check, and canonical ordering so results do not depend on the
// input order or on the standard library's sort (see the tie-break rule in README.md).
#include "detail.hpp"

#include <algorithm>
#include <utility>

namespace lex::detail {

bool sort_values(std::vector<double>& v) {
    // -0.0 + 0.0 is +0.0 and every other value is unchanged, so -0 and +0 (which compare equal)
    // become the same bits and the sorted array is fully determined by the input multiset.
    for (double& x : v) x = x + 0.0;
    if (std::is_sorted(v.begin(), v.end())) return false;
    std::sort(v.begin(), v.end());
    return true;
}

bool sort_pairs(std::vector<double>& x, std::vector<double>& w) {
    const std::size_t n = x.size();
    for (std::size_t i = 0; i < n; ++i) x[i] = x[i] + 0.0;
    auto less = [](double xa, double wa, double xb, double wb) { return xa < xb || (xa == xb && wa < wb); };
    bool sorted = true;
    for (std::size_t i = 1; i < n && sorted; ++i)
        if (less(x[i], w[i], x[i - 1], w[i - 1])) sorted = false;
    if (sorted) return false;
    std::vector<std::pair<double, double>> z(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = {x[i], w[i]};
    // Lexicographic (value, weight): equal values with different weights get one fixed order, so
    // the double-double prefix sums, and hence every cost, are the same on every platform.
    std::sort(z.begin(), z.end(), [&](const auto& a, const auto& b) { return less(a.first, a.second, b.first, b.second); });
    for (std::size_t i = 0; i < n; ++i) { x[i] = z[i].first; w[i] = z[i].second; }
    return true;
}

} // namespace lex::detail
