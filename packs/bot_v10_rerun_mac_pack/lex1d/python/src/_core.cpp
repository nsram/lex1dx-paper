// nanobind bindings: numpy float64 arrays in, a dict of numpy arrays and scalars out.
#include "lex/lex.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>

#include <cstring>
#include <optional>
#include <vector>

namespace nb = nanobind;
using Vec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

namespace {

template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>> to_numpy(const std::vector<T>& v) {
    T* data = new T[v.size() ? v.size() : 1];
    if (!v.empty()) std::memcpy(data, v.data(), v.size() * sizeof(T));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<T*>(p); });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(data, {v.size()}, owner);
}

const char* start_name(lex::Start s) {
    switch (s) {
        case lex::Start::Collapse: return "collapse";
        case lex::Start::EqualFrequency: return "equal_frequency";
        case lex::Start::EqualWidth: return "equal_width";
        case lex::Start::Exact: return "exact";
        case lex::Start::CoarseExact: return "coarse_exact";
    }
    return "?";
}

nb::dict to_dict(const lex::Result& r) {
    nb::dict d;
    std::vector<long long> bounds(r.bounds.begin(), r.bounds.end()), sizes(r.sizes.begin(), r.sizes.end());
    d["bounds"] = to_numpy(bounds);
    d["sizes"] = to_numpy(sizes);
    d["centers"] = to_numpy(r.centers);
    d["upper"] = to_numpy(r.upper);
    d["sse"] = r.sse;
    d["start"] = start_name(r.start);
    d["atoms"] = r.atoms;
    d["sweeps"] = r.sweeps;
    d["moves"] = r.moves;
    d["presorted"] = r.presorted;
    return d;
}

nb::dict lex_fit(Vec x, int k, std::optional<Vec> w, double precision, bool multi_start, double tolerance, std::size_t coarse_groups) {
    lex::Options o;
    o.precision = precision;
    o.multi_start = multi_start;
    o.tolerance = tolerance;
    o.coarse_groups = coarse_groups;
    if (w && w->shape(0) != x.shape(0)) throw std::invalid_argument("weights size must match input");
    lex::Result r;
    {
        nb::gil_scoped_release release;
        r = w ? lex::lex_hybrid(x.data(), x.shape(0), k, w->data(), o) : lex::lex_hybrid(x.data(), x.shape(0), k, o);
    }
    return to_dict(r);
}

nb::dict ckmeans_exact(Vec x, int k, std::optional<Vec> w) {
    if (w && w->shape(0) != x.shape(0)) throw std::invalid_argument("weights size must match input");
    lex::Result r;
    {
        nb::gil_scoped_release release;
        r = w ? lex::ckmeans_exact(x.data(), x.shape(0), k, w->data()) : lex::ckmeans_exact(x.data(), x.shape(0), k);
    }
    return to_dict(r);
}

} // namespace

NB_MODULE(_core, m) {
    m.doc() = "LEX C++ core (internal; use the lex package)";
    m.def("lex", &lex_fit, nb::arg("x"), nb::arg("k"), nb::arg("weights").none(), nb::arg("precision"), nb::arg("multi_start"),
          nb::arg("tolerance"), nb::arg("coarse_groups"));
    m.def("ckmeans_exact", &ckmeans_exact, nb::arg("x"), nb::arg("k"), nb::arg("weights").none());
    lex::Options d;
    m.attr("DEFAULT_COARSE_GROUPS") = d.coarse_groups;
    m.attr("DEFAULT_TOLERANCE") = d.tolerance;
}
