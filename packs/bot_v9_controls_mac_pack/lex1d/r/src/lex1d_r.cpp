// Rcpp binding: one entry point over lex::lex_hybrid / lex_plain / ckmeans_exact. Input checks live in
// R/ckmeans.R.
#include <Rcpp.h>

#include "lex/lex.hpp"

#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

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

} // namespace

// [[Rcpp::export(.lex_fit)]]
Rcpp::List lex_fit(Rcpp::NumericVector x, int k, Rcpp::Nullable<Rcpp::NumericVector> weights, std::string method,
                   double precision, bool multi_start, double tolerance, double coarse_groups) {
    const std::size_t n = x.size();
    const double* w = nullptr;
    Rcpp::NumericVector wv;
    if (weights.isNotNull()) {
        wv = Rcpp::NumericVector(weights);
        if (static_cast<std::size_t>(wv.size()) != n) throw std::invalid_argument("weights size must match input");
        w = wv.begin();
    }
    lex::Options o;
    o.precision = precision;
    o.multi_start = multi_start;
    o.tolerance = tolerance;
    o.coarse_groups = static_cast<std::size_t>(coarse_groups);
    lex::Result r;
    if (method == "ckmeans_exact") r = lex::ckmeans_exact(x.begin(), n, k, w);
    else if (method == "lex_plain") r = lex::lex_plain(x.begin(), n, k, w, o);
    else if (method == "lex_hybrid") r = lex::lex_hybrid(x.begin(), n, k, w, o);
    else throw std::invalid_argument("unknown method " + method);
    // Cluster of each input value (1-based, input order) and per-cluster weighted sums of squares.
    const int kk = r.k();
    Rcpp::IntegerVector cluster(n);
    std::vector<double> withinss(kk, 0.0), wsum(kk, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const int c = r.assign(x[i]);
        cluster[i] = c + 1;
        const double wi = w ? w[i] : 1.0, d = x[i] - r.centers[c];
        withinss[c] += wi * d * d;
        wsum[c] += wi;
    }
    return Rcpp::List::create(
        Rcpp::Named("cluster") = cluster, Rcpp::Named("centers") = r.centers,
        Rcpp::Named("withinss") = withinss, Rcpp::Named("wsum") = wsum,
        Rcpp::Named("count") = r.sizes, Rcpp::Named("upper") = r.upper, Rcpp::Named("sse") = r.sse,
        Rcpp::Named("start") = start_name(r.start), Rcpp::Named("atoms") = r.atoms,
        Rcpp::Named("sweeps") = r.sweeps, Rcpp::Named("moves") = r.moves,
        Rcpp::Named("presorted") = r.presorted);
}

// C strtod (correctly rounded on the supported platforms), for the tests: R's own parser (scan(),
// as.numeric()) can be 1 ulp off where long double is double (Apple Silicon). NA if a string is not
// one whole number.
// [[Rcpp::export(.strtod)]]
Rcpp::NumericVector strtod_(Rcpp::CharacterVector s) {
    Rcpp::NumericVector out(s.size());
    for (R_xlen_t i = 0; i < s.size(); ++i) {
        const char* p = s[i];
        char* end = nullptr;
        const double v = std::strtod(p, &end);
        while (*end == ' ' || *end == '\t' || *end == '\r') ++end;
        out[i] = (end == p || *end != '\0') ? NA_REAL : v;
    }
    return out;
}
