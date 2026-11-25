// Convert one decimal number per line to little-endian IEEE-754 float64.
// Usage: txt_to_f64 IN.txt OUT.f64
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: txt_to_f64 IN.txt OUT.f64\n");
        return 2;
    }
    FILE* in = std::fopen(argv[1], "rb");
    FILE* out = std::fopen(argv[2], "wb");
    if (!in || !out) {
        std::fprintf(stderr, "cannot open input or output\n");
        return 2;
    }
    std::string line;
    line.reserve(64);
    int c;
    std::size_t n = 0;
    while (c = std::fgetc(in), c != EOF) {
        if (c == '\n') {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            const std::size_t s = line.find_first_not_of(" \t");
            if (s != std::string::npos && line[s] != '#') {
                char* end = nullptr;
                const double v = std::strtod(line.c_str() + s, &end);
                if (end == line.c_str() + s) {
                    std::fprintf(stderr, "not a number at line after %zu values: %s\n", n, line.c_str());
                    return 2;
                }
                if (std::fwrite(&v, 8, 1, out) != 1) return 2;
                ++n;
            }
            line.clear();
        } else {
            line.push_back(static_cast<char>(c));
        }
    }
    if (!line.empty()) {
        if (line.back() == '\r') line.pop_back();
        const std::size_t s = line.find_first_not_of(" \t");
        if (s != std::string::npos && line[s] != '#') {
            const double v = std::strtod(line.c_str() + s, nullptr);
            if (std::fwrite(&v, 8, 1, out) != 1) return 2;
            ++n;
        }
    }
    std::fclose(in);
    std::fclose(out);
    std::printf("WROTE n=%zu path=%s\n", n, argv[2]);
    return n == 0 ? 2 : 0;
}
