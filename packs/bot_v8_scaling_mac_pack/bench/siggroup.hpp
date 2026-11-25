// 5-significant-digit grouping of an already-sorted array.
// Reference key is the run-5 / precgroup nearbyint rule.
// The fast path is a linear neighbour scan: within a conservative decade it
// reuses the scale and only calls nearbyint; it never hashes and does not
// allocate a per-element key array. Groups are identical to the reference
// (same boundaries, therefore same count, sum, and sumsq).
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace siggroup {

inline double sig_key(double x, int s) {
    if (x == 0.0) return 0.0;
    const int e = static_cast<int>(std::floor(std::log10(std::fabs(x))));
    const double q = std::pow(10.0, static_cast<double>(e - s + 1));
    return std::nearbyint(x / q) * q;
}

struct Groups {
    std::vector<int> start;     // length D'+1, last = n
    std::vector<double> mean, w;
    std::vector<long double> sum, sumsq;
};

inline Groups group_ref(const std::vector<double>& sorted, int sig) {
    Groups g;
    const std::size_t n = sorted.size();
    g.start.reserve(n / 8 + 2);
    g.mean.reserve(n / 8 + 2);
    g.w.reserve(n / 8 + 2);
    g.sum.reserve(n / 8 + 2);
    g.sumsq.reserve(n / 8 + 2);
    std::size_t i = 0;
    while (i < n) {
        const double kk = sig_key(sorted[i], sig);
        std::size_t j = i;
        long double s = 0, q = 0;
        while (j < n && sig_key(sorted[j], sig) == kk) {
            const long double x = sorted[j];
            s += x;
            q += x * x;
            ++j;
        }
        g.start.push_back(static_cast<int>(i));
        g.sum.push_back(s);
        g.sumsq.push_back(q);
        g.mean.push_back(static_cast<double>(s / (j - i)));
        g.w.push_back(static_cast<double>(j - i));
        i = j;
    }
    g.start.push_back(static_cast<int>(n));
    return g;
}

// Linear pass over a sorted array. For each group we compute the 5-sig key
// once, then consume the run of neighbours that share it:
//   * identical values (no key work),
//   * conservative interior of the same decimal decade: compare x against a
//     slightly shrunk upper bound of the nearbyint slot (no log10/pow),
//   * decade edges, zeros, and sign changes: fall back to sig_key.
// No hash table and no per-element key array. Groups match group_ref.
inline Groups group_fast(const std::vector<double>& sorted, int sig) {
    Groups g;
    const std::size_t n = sorted.size();
    g.start.reserve(n / 8 + 2);
    g.mean.reserve(n / 8 + 2);
    g.w.reserve(n / 8 + 2);
    g.sum.reserve(n / 8 + 2);
    g.sumsq.reserve(n / 8 + 2);

    const double ninf = -__builtin_inf();
    std::size_t i = 0;
    while (i < n) {
        const double xi = sorted[i];
        const double kk = sig_key(xi, sig);
        const bool zero = (xi == 0.0);
        double mag_lo = 0, mag_hi = 0, hi = 0, qscale = 0, slot = 0;
        if (!zero) {
            const int e = static_cast<int>(std::floor(std::log10(std::fabs(xi))));
            qscale = std::pow(10.0, static_cast<double>(e - sig + 1));
            slot = std::nearbyint(xi / qscale);
            const double p_lo = std::pow(10.0, static_cast<double>(e));
            const double p_hi = std::pow(10.0, static_cast<double>(e + 1));
            mag_lo = p_lo * 1.0000000001;
            mag_hi = p_hi * 0.9999999999;
            if (mag_lo >= mag_hi) { mag_lo = p_hi; mag_hi = p_lo; } // empty fast range
            hi = (slot + 0.5) * qscale;
            for (int u = 0; u < 8; ++u) hi = std::nextafter(hi, ninf);
        }
        auto same_as = [&](double xj) -> bool {
            if (xj == xi) return true;
            if (zero || xj == 0.0 || ((xj > 0) != (xi > 0))) return sig_key(xj, sig) == kk;
            const double ax = std::fabs(xj);
            if (ax >= mag_lo && ax <= mag_hi) {
                if (xj > 0 && xj <= hi) return true;
                return std::nearbyint(xj / qscale) == slot;
            }
            return sig_key(xj, sig) == kk;
        };
        std::size_t j = i;
        long double s = 0, qq = 0;
        while (j < n && same_as(sorted[j])) {
            const long double x = sorted[j];
            s += x;
            qq += x * x;
            ++j;
        }
        g.start.push_back(static_cast<int>(i));
        g.sum.push_back(s);
        g.sumsq.push_back(qq);
        g.mean.push_back(static_cast<double>(s / (j - i)));
        g.w.push_back(static_cast<double>(j - i));
        i = j;
    }
    g.start.push_back(static_cast<int>(n));
    return g;
}

inline bool identical(const Groups& a, const Groups& b) {
    if (a.start != b.start) return false;
    if (a.w != b.w) return false;
    if (a.sum.size() != b.sum.size() || a.sumsq.size() != b.sumsq.size()) return false;
    for (std::size_t i = 0; i < a.sum.size(); ++i) {
        if (a.sum[i] != b.sum[i]) return false;
        if (a.sumsq[i] != b.sumsq[i]) return false;
    }
    // means are derived from the same long-double sum and count
    if (a.mean.size() != b.mean.size()) return false;
    for (std::size_t i = 0; i < a.mean.size(); ++i) {
        if (a.mean[i] != b.mean[i]) return false;
    }
    return true;
}

} // namespace siggroup
