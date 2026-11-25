// 12 values, k = 3: same output format as examples/twelve.py and the Scala examples.Twelve.
#include <lex/lex.hpp>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

static unsigned long long bits(double v) { std::uint64_t u; std::memcpy(&u, &v, sizeof u); return u; }

int main() {
    const std::vector<double> x{5.5, 1.2, 9.8, 2.0, 1.0, 5.0, 10.4, 1.9, 6.1, 1.2, 2.1, 5.5};
    const lex::Result r = lex::lex_hybrid(x, 3);
    std::printf("sizes  ");
    for (int s : r.sizes) std::printf(" %d", s);
    std::printf("\ncenters");
    for (double c : r.centers) std::printf(" %.6f", c);
    std::printf("\nsse     %.6f  bits %016llx\ncenters bits", r.sse, bits(r.sse));
    for (double c : r.centers) std::printf(" %016llx", bits(c));
    std::printf("\n");
    return 0;
}
