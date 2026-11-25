// Tiny driver for gronlund/kmeans1d. Reads one double per line, plus k.
// Sorts a copy (the library assumes sorted input). Prints SSE, centers, path, time.
//
//   gronlund_driver --alg ALG --k K [--mode compute|report] [--out FILE] DATAFILE
//
// ALG: linear, monotone, slow, hirc, wilber-binary, wilber-interp
// mode compute: cost only for DP methods (two-row fill). report: compute_and_report
// (Hirschberg-style reconstruction for linear/monotone; same as compute for wilber/hirc).

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "kmeans.h"

static void usage() {
    std::fprintf(stderr,
        "usage: gronlund_driver --alg ALG --k K [--mode compute|report] [--out FILE] DATAFILE\n");
    std::exit(2);
}

static std::vector<double> read_points(const char* path) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot read %s\n", path);
        std::exit(2);
    }
    std::vector<double> xs;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto s = line.find_first_not_of(" \t");
        if (s == std::string::npos || line[s] == '#') continue;
        char* end = nullptr;
        const char* p = line.c_str() + s;
        double v = std::strtod(p, &end);
        if (end == p) continue;
        xs.push_back(v);
    }
    return xs;
}

int main(int argc, char** argv) {
    std::string alg, mode = "report", out_path, data_path;
    std::size_t k = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--alg" && i + 1 < argc) alg = argv[++i];
        else if (a == "--k" && i + 1 < argc) k = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        else if (a == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (a == "--out" && i + 1 < argc) out_path = argv[++i];
        else if (a[0] != '-') data_path = a;
        else usage();
    }
    if (alg.empty() || k == 0 || data_path.empty()) usage();
    if (mode != "compute" && mode != "report") usage();

    std::vector<double> points = read_points(data_path.c_str());
    const int already_sorted = std::is_sorted(points.begin(), points.end()) ? 1 : 0;
    if (!already_sorted) std::sort(points.begin(), points.end());
    const std::size_t n = points.size();
    if (n == 0) {
        std::fprintf(stderr, "empty input\n");
        return 2;
    }

    std::unique_ptr<kmeans> fn;
    if (alg == "linear") fn.reset(new kmeans_linear(points));
    else if (alg == "monotone") fn.reset(new kmeans_monotone(points));
    else if (alg == "slow") fn.reset(new kmeans_slow(points));
    else if (alg == "hirc") fn.reset(new kmeans_hirschberg_larmore(points));
    else if (alg == "wilber-binary") {
        auto w = new kmeans_wilber(points);
        w->set_search_strategy(search_strategy::BINARY);
        fn.reset(w);
    } else if (alg == "wilber-interp") {
        auto w = new kmeans_wilber(points);
        w->set_search_strategy(search_strategy::INTERPOLATION);
        fn.reset(w);
    } else {
        std::fprintf(stderr, "unknown alg: %s\n", alg.c_str());
        return 2;
    }

    const auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<kmeans_result> res =
        (mode == "report") ? fn->compute_and_report(k) : fn->compute(k);
    const auto t1 = std::chrono::steady_clock::now();
    const double time_s = std::chrono::duration<double>(t1 - t0).count();

    std::ostringstream oss;
    oss << std::setprecision(17);
    oss << "RESULT"
        << " alg=" << alg
        << " n=" << n
        << " k=" << k
        << " mode=" << mode
        << " sse=" << (res ? res->cost : -1)
        << " n_centers=" << (res ? res->centers.size() : 0)
        << " n_path=" << (res ? res->path.size() : 0)
        << " time_s=" << std::setprecision(9) << time_s
        << " sorted_input=" << already_sorted
        << " lib_name=" << fn->name()
        << "\n";
    oss << "CENTERS";
    if (res) {
        for (double c : res->centers) oss << " " << std::setprecision(17) << c;
    }
    oss << "\nPATH";
    if (res) {
        for (std::size_t p : res->path) oss << " " << p;
    }
    oss << "\n";
    const std::string text = oss.str();
    std::fwrite(text.data(), 1, text.size(), stdout);
    std::fflush(stdout);
    if (!out_path.empty()) {
        std::ofstream out(out_path);
        out << text;
    }
    return 0;
}
