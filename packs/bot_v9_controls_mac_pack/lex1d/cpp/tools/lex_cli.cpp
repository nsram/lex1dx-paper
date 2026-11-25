// lex_cli: fit one column of numbers and print the result in a stable, bit-exact text form.
//
//   lex_cli FILE K [--preset hybrid|plain|deck|exact] [--coarse M] [--weights WFILE]
//
// Presets: hybrid = lex::lex_hybrid (the default; --coarse M sets its coarse_groups), plain =
// lex::lex_plain, deck = lex::lex_plain with lex::deck(), exact = lex::ckmeans_exact.
//
// Input (FILE, and WFILE with one weight per line):
//   - one "value" or "value,weight" per line (spaces around the comma are fine); every data line
//     of a file must have the same form; "value,weight" cannot be combined with --weights;
//   - blank lines and lines whose first non-space character is '#' are skipped; "\r\n" line
//     endings are accepted;
//   - the first remaining line is skipped as a header if it does not start with a number;
//   - numbers are parsed with strtod (decimal or hex floats; "nan" and "inf" parse but are then
//     rejected by the library with error=...); anything else on a line is an error (exit status 2).
// Output, one "key=value" per line: bounds, sizes, sse (%.17g), sse_bits and centers_bits
// (IEEE-754 bit patterns in hex), start, atoms, sweeps, moves. Used by the Python and Scala
// cross-checks; fitting errors print "error=<message>" and exit with status 1.
#include "lex/lex.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

[[noreturn]] void bad_input(const char* path, int line, const std::string& msg) {
    std::fprintf(stderr, "%s:%d: %s\n", path, line, msg.c_str());
    std::exit(2);
}

/// Parses a number starting at p (after optional spaces); returns false if there is none.
bool number(const char*& p, double& v) {
    while (*p == ' ' || *p == '\t') ++p;
    char* end = nullptr;
    v = std::strtod(p, &end);
    if (end == p) return false;
    p = end;
    while (*p == ' ' || *p == '\t') ++p;
    return true;
}

/// Reads values (and weights, if the lines are "value,weight") from a text file; see the header.
void read_input(const char* path, std::vector<double>& xs, std::vector<double>& ws, bool allow_pairs) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(2); }
    std::string line;
    int lineno = 0, form = 0; // form: 0 unknown, 1 "value", 2 "value,weight"
    bool first = true;
    while (std::getline(in, line)) {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::size_t s = line.find_first_not_of(" \t");
        if (s == std::string::npos || line[s] == '#') continue;
        const char* p = line.c_str() + s;
        double v = 0, w = 0;
        if (!number(p, v)) {
            if (first) { first = false; continue; } // header
            bad_input(path, lineno, "not a number: " + line);
        }
        first = false;
        int f = 1;
        if (*p == ',') {
            ++p;
            if (!allow_pairs) bad_input(path, lineno, "expected one number per line: " + line);
            if (!number(p, w)) bad_input(path, lineno, "expected value,weight: " + line);
            f = 2;
        }
        if (*p != '\0') bad_input(path, lineno, "unexpected text after the number(s): " + line);
        if (form != 0 && f != form) bad_input(path, lineno, "mixes \"value\" and \"value,weight\" lines");
        form = f;
        xs.push_back(v);
        if (f == 2) ws.push_back(w);
    }
}

std::uint64_t bits(double x) { std::uint64_t u; std::memcpy(&u, &x, sizeof u); return u; }

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

int usage() {
    std::fprintf(stderr, "usage: lex_cli FILE K [--preset hybrid|plain|deck|exact] [--coarse M] [--weights WFILE]\n");
    return 2;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) return usage();
    std::vector<double> xs, ws, none;
    read_input(argv[1], xs, ws, true);
    bool weighted = !ws.empty();
    const int k = std::atoi(argv[2]);
    std::string preset = "hybrid";
    long coarse = -1;
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (i + 1 >= argc) return usage();
        if (a == "--preset") preset = argv[++i];
        else if (a == "--coarse") coarse = std::atol(argv[++i]);
        else if (a == "--weights") {
            if (weighted) { std::fprintf(stderr, "--weights given but FILE already has value,weight lines\n"); return 2; }
            read_input(argv[++i], ws, none, false);
            weighted = true;
        }
        else return usage();
    }
    if (preset != "hybrid" && preset != "plain" && preset != "deck" && preset != "exact") return usage();
    lex::Options o = preset == "deck" ? lex::deck() : lex::Options{};
    if (coarse >= 0) o.coarse_groups = static_cast<std::size_t>(coarse);
    const double* w = nullptr;
    lex::Result r;
    try {
        if (weighted && ws.size() != xs.size()) throw std::invalid_argument("weights size must match input");
        if (weighted) w = ws.data();
        if (preset == "exact") r = lex::ckmeans_exact(xs.data(), xs.size(), k, w);
        else if (preset == "hybrid") r = lex::lex_hybrid(xs.data(), xs.size(), k, w, o);
        else r = lex::lex_plain(xs.data(), xs.size(), k, w, o);
    } catch (const std::exception& e) {
        std::printf("error=%s\n", e.what());
        return 1;
    }
    auto ints = [](const std::vector<int>& v) {
        std::string s;
        for (std::size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + std::to_string(v[i]);
        return s;
    };
    std::printf("bounds=%s\nsizes=%s\nsse=%.17g\nsse_bits=%016llx\ncenters_bits=", ints(r.bounds).c_str(), ints(r.sizes).c_str(),
                r.sse, static_cast<unsigned long long>(bits(r.sse)));
    for (std::size_t i = 0; i < r.centers.size(); ++i)
        std::printf("%s%016llx", i ? "," : "", static_cast<unsigned long long>(bits(r.centers[i])));
    std::printf("\nstart=%s\natoms=%d\nsweeps=%d\nmoves=%d\n", start_name(r.start), r.atoms, r.sweeps, r.moves);
    return 0;
}
