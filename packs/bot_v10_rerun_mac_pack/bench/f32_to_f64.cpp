// Promote little-endian float32 values to little-endian float64 (exact).
//   f32_to_f64 IN.f32le OUT.f64 [N]
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: f32_to_f64 IN.f32le OUT.f64 [N]\n");
        return 2;
    }
    FILE* in = std::fopen(argv[1], "rb");
    if (!in) { std::fprintf(stderr, "cannot read %s\n", argv[1]); return 2; }
    std::fseek(in, 0, SEEK_END);
    const long bytes = std::ftell(in);
    std::fseek(in, 0, SEEK_SET);
    if (bytes < 0 || bytes % 4 != 0) { std::fprintf(stderr, "bad f32 file %s\n", argv[1]); return 2; }
    std::size_t n = static_cast<std::size_t>(bytes) / 4;
    if (argc > 3) {
        const std::size_t cap = static_cast<std::size_t>(std::strtoull(argv[3], nullptr, 10));
        if (cap < n) n = cap;
    }
    std::vector<float> x(n);
    if (std::fread(x.data(), 4, n, in) != n) { std::fprintf(stderr, "short read %s\n", argv[1]); return 2; }
    std::fclose(in);
    std::vector<double> y(n);
    for (std::size_t i = 0; i < n; ++i) y[i] = static_cast<double>(x[i]);
    FILE* out = std::fopen(argv[2], "wb");
    if (!out) { std::fprintf(stderr, "cannot write %s\n", argv[2]); return 2; }
    if (std::fwrite(y.data(), 8, n, out) != n) { std::fprintf(stderr, "short write %s\n", argv[2]); return 2; }
    std::fclose(out);
    std::printf("WROTE path=%s n=%zu\n", argv[2], n);
    return 0;
}
