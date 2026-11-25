#include "detail.hpp"

namespace lex::detail {

Atoms::Atoms(std::vector<double> sorted_values, double precision) : sorted(std::move(sorted_values)) {
    n = static_cast<int>(sorted.size());
    double s = 0;
    for (int i = 0; i < n; ++i) s += sorted[i];
    const double m = s / n;
    double c = 0;
    for (int i = 0; i < n; ++i) c += sorted[i] - m;
    center = m + c / n;

    atom_id.resize(n);
    start.resize(n + 1);
    value.resize(n);
    int i = 0;
    while (i < n) {
        const int first = i;
        double sum = sorted[i] - center;
        atom_id[i] = d;
        ++i;
        // Equal values always share an atom (tie rule), whatever the precision test says.
        while (i < n && (sorted[i] == sorted[i - 1] ||
                         (precision > 0 && sorted[i] - center - sum / (i - first) < precision))) {
            sum += sorted[i] - center;
            atom_id[i] = d;
            ++i;
        }
        start[d] = first;
        value[d] = precision == 0 ? sorted[first] - center : sum / (i - first);
        ++d;
    }
    start[d] = n;
    start.resize(d + 1);
    value.resize(d);

    prefix_hi.assign(n + 1, 0.0);
    prefix_lo.assign(n + 1, 0.0);
    double hi = 0, lo = 0;
    for (int e = 0; e < n; ++e) {
        const double x = value[atom_id[e]];
        const double t = hi + x; // TwoSum: t + err == hi + x exactly
        const double bv = t - hi;
        lo += (hi - (t - bv)) + (x - bv);
        hi = t;
        prefix_hi[e + 1] = hi;
        prefix_lo[e + 1] = lo;
    }
}

double Atoms::sse(const std::vector<int>& b) const {
    if (weighted) {
        double total = 0;
        for (std::size_t i = 0; i + 1 < b.size(); ++i) total += range_sse(b[i], b[i + 1]);
        return total;
    }
    double total = 0;
    for (std::size_t i = 0; i + 1 < b.size(); ++i) {
        const double m = mean(b[i], b[i + 1]);
        for (int e = b[i]; e < b[i + 1]; ++e) {
            const double dv = value[atom_id[e]] - m;
            total += dv * dv;
        }
    }
    return total;
}

namespace {
void two_sum(std::vector<double>& hi, std::vector<double>& lo, int j, double x) {
    const double s = hi[j] + x;
    const double bv = s - hi[j];
    hi[j + 1] = s;
    lo[j + 1] = lo[j] + (hi[j] - (s - bv)) + (x - bv);
}
} // namespace

Atoms::Atoms(std::vector<double> sorted_values, std::vector<double> sorted_weights, double precision)
    : sorted(std::move(sorted_values)), weighted(true), elem_weight(std::move(sorted_weights)) {
    n = static_cast<int>(sorted.size());
    double W = 0, S = 0;
    for (int i = 0; i < n; ++i) { W += elem_weight[i]; S += elem_weight[i] * sorted[i]; }
    const double m = S / W;
    double c = 0;
    for (int i = 0; i < n; ++i) c += elem_weight[i] * (sorted[i] - m);
    center = m + c / W;

    atom_id.resize(n);
    start.resize(n + 1);
    value.resize(n);
    int i = 0;
    while (i < n) {
        const int first = i;
        double wrun = elem_weight[i], srun = elem_weight[i] * (sorted[i] - center);
        int crun = 1;
        double s_unw = sorted[i] - center;
        atom_id[i] = d;
        ++i;
        while (i < n && (sorted[i] == sorted[i - 1] ||
                         (precision > 0 && (wrun > 0 ? sorted[i] - center - srun / wrun
                                                     : sorted[i] - center - s_unw / crun) < precision))) {
            wrun += elem_weight[i];
            srun += elem_weight[i] * (sorted[i] - center);
            s_unw += sorted[i] - center;
            ++crun;
            atom_id[i] = d;
            ++i;
        }
        start[d] = first;
        if (precision == 0) value[d] = sorted[first] - center;
        else value[d] = wrun > 0 ? srun / wrun : s_unw / crun;
        ++d;
    }
    start[d] = n;
    start.resize(d + 1);
    value.resize(d);

    w_hi.assign(n + 1, 0.0); w_lo.assign(n + 1, 0.0);
    ws_hi.assign(n + 1, 0.0); ws_lo.assign(n + 1, 0.0);
    wq_hi.assign(n + 1, 0.0); wq_lo.assign(n + 1, 0.0);
    for (int e = 0; e < n; ++e) {
        const double w = elem_weight[e];
        const double z = value[atom_id[e]];
        two_sum(w_hi, w_lo, e, w);
        two_sum(ws_hi, ws_lo, e, w * z);
        two_sum(wq_hi, wq_lo, e, w * z * z);
    }
}

} // namespace lex::detail
