# lex1d

Repository: https://github.com/nsram/lex1d


LEX partitions one-dimensional data into k contiguous bins and keeps the within-bin sum of squared
deviations low. This is 1-D k-means, also called optimal scalar quantization. Typical uses are
quantizing model weights and embeddings, building histograms and Jenks-style class breaks, and
discretizing a column.

LEX does this in three steps:

1. **Sort the data and group equal values into atoms.** An atom is never split.
2. **Collapse:** greedily merge the cheapest adjacent pair of bins until k bins remain.
3. **Refine:** improve the bins with 3-bin-window moves until no move helps.

Every language offers the same three methods:

| Method | What it does |
|---|---|
| `lex_hybrid` | The default. An exact dynamic program on a coarse 4,096-group summary seeds Refine (see "Accuracy and speed" below). |
| `lex_plain` | Plain LEX: Collapse seeds Refine. |
| `ckmeans_exact` | The exact dynamic program on all the data (the Ckmeans.1d.dp solver), no LEX steps. |

The 2024 deck preset (a 1e-4 grouping grid and three starts) is `lex_plain` with `precision = 1e-4`
and `multi_start` on.

The same algorithm is provided in four forms:

| Directory | What it holds |
|---|---|
| `cpp/` | The C++17 library: CMake, ctest and a small `lex_cli` tool |
| `python/` | Python bindings (nanobind), using numpy arrays in and out, with optional weights |
| `scala/` | The Scala 3 reference implementation, with weights, cross-checked bit-for-bit against C++ |
| `r/` | R package `lex1d` (Rcpp) on the C++ core; `ckmeans()` (= `lex_hybrid()`) follows the interface of Ckmeans.1d.dp |

The math is in [`docs/algorithm.md`](docs/algorithm.md).

## Install

Exact, tested build and test steps for Linux, macOS and Windows are in [BUILDING.md](BUILDING.md).
In short, from a checkout of the repository:

```sh
git clone https://github.com/nsram/lex1d && cd lex1d
python fixtures/gen_fixtures.py                     # check the committed test fixtures
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release && ctest --test-dir build -C Release
pip install ./python && pytest python/tests         # Python package lex1d, import name lex
cd scala && sbt test && cd ..                       # Scala (sbt, JDK 17)
R CMD INSTALL r                                     # R package lex1d (needs Rcpp); tests: see r/README.md
```

To use the C++ library from your own CMake project: `add_subdirectory(lex1d/cpp)` and link
`lex::lex` (see `examples/cpp`). The Python distribution is named `lex1d` (like the repository) and
the import name is `lex`. It is not on PyPI yet. The R package is also named `lex1d` and is not on
CRAN yet; see [r/README.md](r/README.md).

## Quick start

C++:

```cpp
#include <lex/lex.hpp>
lex::Result r = lex::lex_hybrid(values, 16);              // the default: coarse exact start
lex::Result p = lex::lex_plain(values, 16);               // plain LEX (fastest at mid N)
lex::Result d = lex::lex_plain(values, 16, lex::deck());  // 2024 deck preset: 1e-4 grid
lex::Result e = lex::ckmeans_exact(values, 16);           // global optimum (DP)
lex::Result w = lex::lex_hybrid(values, weights, 16);     // weighted
// r.bounds (k + 1 cuts into the sorted data), r.sizes, r.centers, r.upper, r.sse, r.assign(x)
```

Python:

```python
import numpy as np, lex
x = np.random.default_rng(0).lognormal(size=100_000)
r = lex.lex_hybrid(x, 16)                   # the default
p = lex.lex_plain(x, 16)                    # deck: lex.lex_plain(x, 16, precision=1e-4, multi_start=True)
e = lex.ckmeans_exact(x, 16)
w = lex.lex_hybrid(x, 16, weights=np.ones_like(x))
labels = r.assign(x)                        # bin index per value (numpy array)
print(r.sse, e.sse, r.centers)
```

Scala:

```scala
import lex.{Lex, Options}
val r = Lex.lexHybrid(xs, 16)               // the default
val p = Lex.lexPlain(xs, 16)                // deck: Lex.lexPlain(xs, 16, Options.deck)
val e = Lex.ckmeansExact(xs, 16)
val w = Lex.lexHybrid(xs, weights, 16)      // weighted; Lex.ckmeansExact(xs, weights, 16)
println((r.sse, e.sse, r.counts.toSeq, r.assign(0.5)))
```

R:

```r
library(lex1d)
x <- rlnorm(100000)
r <- lex_hybrid(x, 16)                      # the default; ckmeans(x, 16) is the same
p <- lex_plain(x, 16)                       # deck: lex_plain(x, 16, precision = 1e-4, multi_start = TRUE)
e <- ckmeans_exact(x, 16)
w <- ckmeans(x, 16, y = rep(1, length(x)))  # weighted (y, as in Ckmeans.1d.dp)
r$cluster; r$centers; r$size; r$tot.withinss   # Ckmeans.1d.dp-style result
```

## Same algorithm, three languages

LEX sorts the values, groups equal values into atoms, and seeds the bins from an exact dynamic
program over a coarse summary of at most 4,096 groups (or greedy merging with `lex_plain`). Then it
moves single atoms across bin boundaries, three bins at a time, while that lowers the sum of squared
errors. Ties between equal-SSE answers go to the largest cuts, as in Ckmeans.1d.dp. The details are
in [docs/algorithm.md](docs/algorithm.md).

C++, Python, R and Scala implement this with the same arithmetic (Python and R call the C++ core). They return bit-identical cuts,
SSE and centres on all 356 fixture runs (89 cases × 4 settings: hybrid, plain, deck, exact; including 64 weighted and 16 with k
clamped). Three tests check this: `python/tests/test_lex.py::test_bit_identical_to_cpp_cli` (Python vs the
C++ CLI) and `scala/src/test/scala/lex/CrossCheckTest.scala` (Scala vs the C++ output) and `r/tests/testthat/test-fixtures.R` (R vs the C++ output).

Example: 12 values, k = 3 (`examples/twelve.py`, `examples/cpp/twelve.cpp`,
`sbt "runMain examples.Twelve"`). Each of the three programs printed exactly this:

```
sizes   6 4 2
centers 1.566667 5.525000 10.100000
sse     1.960833  bits 3fff5f92c5f92c5d
centers bits 3ff9111111111112 401619999999999a 4024333333333334
```

The only differences between the implementations:

| | C++ / Python | Scala |
|---|---|---|
| Presorted input | detected; the sort is skipped (`presorted` in the result) | always sorts; no flag |
| Result field names | `sizes`, `centers` | `counts`, `means` |
| Numbers in error messages | C formatting (`1e-04`) | Java formatting (`1.0E-4`) |

Python (`UserWarning`) and R (`warning()`) warn when they clamp k; C++ and Scala clamp silently. The R
result uses Ckmeans.1d.dp's names (`cluster`, `centers`, `size`, `withinss`, ...); its other differences
from Ckmeans.1d.dp are listed in [r/README.md](r/README.md).

## Input format

| Interface | Values | Weights (optional) |
|---|---|---|
| Python `lex.lex_hybrid(x, k, weights=...)` (and the other two) | any 1-D array-like; converted to contiguous float64 (ints, float32, lists are fine); 2-D input is a `ValueError` | same length, finite, ≥ 0, positive total |
| C++ `lex::lex_hybrid`, `lex_plain`, `ckmeans_exact` | `std::vector<double>`, or `const double*` plus length | `std::vector<double>` of the same length, or a pointer (`nullptr` = unweighted) |
| Scala `Lex.lexHybrid`, `lexPlain`, `ckmeansExact` | `Array[Double]` | `Array[Double]` of the same length |
| R `lex_hybrid(x, k, y = ...)` (and `lex_plain`, `ckmeans_exact`, `ckmeans`) | numeric vector (integer is fine); `NA` is an error | `y`, same length, finite, ≥ 0, positive total |
| `lex_cli FILE K` | text file, see below | `value,weight` lines in FILE, or `--weights WFILE` |

All interfaces:
- reject empty input, NaN and ±inf, and k < 1 (Python `ValueError`, C++ `std::invalid_argument`,
  Scala `IllegalArgumentException`, R `stop()`);
- treat −0.0 as +0.0;
- never modify the input;
- clamp k to the number of distinct values (or deck groups).

`lex_cli` reads one `value` or one `value,weight` per line, with optional spaces around the comma.
Every data line of a file must have the same form, and `value,weight` cannot be combined with
`--weights`. Other rules:
- Blank lines and lines starting with `#` are skipped, and `\r\n` line endings are accepted.
- The first remaining line is skipped as a header if it does not start with a number.
- Numbers are parsed with `strtod`. Any other text on a line stops with `FILE:LINE: message` and
  exit status 2. A fitting error (for example a NaN value) prints `error=...` and exits with
  status 1.

## Examples

Every example generates its own data and downloads nothing. The outputs below are real runs on the
development box (Linux x86-64 VM, single thread); times there are indicative only.

`python examples/quickstart.py` fits uniform, normal and lognormal samples (N = 200,000, 10 % of
values rounded to 2 decimals, so there are ties) with k = 16. Columns: the gap to `ckmeans_exact` for
`lex_hybrid` and `lex_plain`, how many times larger a uniform min-max grid's MSE is, and the times.

```
unweighted
dist              N   k   gap hybrid  gap plain  uniform/LEX   t hybrid  t plain  t exact
uniform      200000  16    0.00e+00%     0.384%        1.00x     77.8ms   75.7ms  130.0ms
normal       200000  16    3.08e-07%     0.018%        2.90x     70.1ms  100.9ms  131.4ms
lognormal    200000  16    1.19e-07%     0.060%       97.64x     91.7ms   81.7ms  130.7ms

weighted (lognormal weights)
dist              N   k   gap hybrid  gap plain  uniform/LEX   t hybrid  t plain  t exact
uniform      200000  16    2.14e-07%     0.016%        1.00x    128.9ms   91.4ms  185.1ms
normal       200000  16    2.62e-06%     0.018%        2.93x     94.8ms   83.9ms  166.6ms
lognormal    200000  16    0.00e+00%     0.596%       63.92x     94.7ms  100.2ms  175.8ms

normal sample, k = 16
cuts (upper value of each bin): [-2.398 -1.836 -1.425 -1.087 -0.791 -0.517 -0.257 -0.003  0.254  0.517  0.794  1.095  1.434  1.844
  2.406  4.523]
centres: [-2.731 -2.065 -1.608 -1.242 -0.933 -0.649 -0.384 -0.13   0.124  0.383  0.651  0.937  1.252  1.615
  2.073  2.74 ]
sizes: [ 1662  5007  8628 12244 15134 17797 19074 20171 20302 19548 17656 15398 12243  8623  4917  1596]  start: coarse_exact  atoms: 180513
```

`python examples/embeddings.py` quantizes each row of a synthetic 2,000 × 64 matrix to 4 bits. The
rows are uniform, normal or heavy-tailed. The table shows MSE per row family:

```
2000 x 64 matrix, per-row k = 16 (4 bits); MSE per row family
method                 uniform        normal  heavy-tailed         all  vs exact     time
lex_hybrid           6.419e-04     7.542e-04     1.904e-03   1.085e-03    1.000x    0.07s
lex_plain            6.531e-04     7.663e-04     1.937e-03   1.104e-03    1.017x    0.04s
ckmeans_exact        6.419e-04     7.542e-04     1.904e-03   1.085e-03    1.000x    0.06s
uniform min-max      1.305e-03     1.934e-03     1.829e-02   6.980e-03    6.432x    0.00s
```

The C++ example (`cmake -S examples/cpp -B build-examples && cmake --build build-examples --config Release && build-examples/quickstart`) and the Scala example (`cd scala && sbt "runMain examples.Quickstart"`) use
the same portable random generator. Both printed the same gaps and ratios; Scala's times include
JIT warm-up. C++ output:

```
dist              N   k   gap hybrid  gap plain  uniform/LEX   t hybrid  t plain  t exact
uniform      200000  16    0.00e+00%     1.303%        1.00x     87.7ms   78.0ms  142.3ms
normal       200000  16    0.00e+00%     0.167%        2.70x     74.3ms   71.9ms  142.2ms
lognormal    200000  16    0.00e+00%     3.168%       32.23x     76.9ms   91.6ms  137.2ms
weighted lognormal: gap hybrid 2.07e-06%, k used 16
```

## Accuracy and speed

Measured on a MacBook Pro with an Apple M3 Max and 128 GB of memory. LEX times are for the C++
library with the default M = 4,096; the baseline is the exact-DP R package Ckmeans.1d.dp. The
benchmark harness and raw results live outside this repository. `bench/plot_timing.py` redraws a
timing figure from a run's `results.csv`. "Excess error" is the relative excess SSE of
`lex_hybrid` over the exact optimum.

- **Exact when D ≤ M.** If the data have at most M distinct values (D ≤ 4,096 by default), the
  coarse start is the exact DP on all the data, and `lex_hybrid` returns `ckmeans_exact`'s
  partition.
- **Accuracy.** On continuous simulated data with up to 1 million points, the excess error grew
  roughly as k²: below 0.002 % for k ≤ 64 and below 0.08 % for k ≤ 512, with no trend in N.
- **Speed.** `lex_hybrid` was faster than Ckmeans.1d.dp in every continuous cell from k = 8, at
  least 2.2× faster for k ≥ 16, and a median 68× faster at N = 1 million, k = 512.
- **Memory.** The input plus a k × M table; the table is at most 8.4 MB for k ≤ 512.
- **Audio, 8.2 million samples, k = 512.** `lex_hybrid` took 3.0 s and 0.61 GB of peak memory with
  0.028 % excess error. Ckmeans.1d.dp took 172 s and 68 GB.

**Choosing M.** The default is `coarse_groups = 4096`. A larger M, for example 8,192, is an option,
not the default: at N = 1 million on continuous data, each doubling of M lowered the excess error
about 3.9× at little time cost. The k × M table grows in proportion to M. An M at or above the
number of distinct values gives the exact optimum.

**Ties.** Among equal-SSE partitions, every solver returns the one with the largest cuts, as
Ckmeans.1d.dp does.

**Reproducibility.** The C++ code is compiled without floating-point contraction (`-ffp-contract=off`
on GCC/Clang, `/fp:precise` on MSVC 2022; see `cpp/cmake/LexFloatingPoint.cmake`). Results are bit-identical
between the C++ CLI and the Python bindings. On the fixture set they are also bit-identical to the
Scala implementation and the R package, weighted cases included. This holds for the same input
doubles: R's own text parser (`scan()`, `read.csv()`) can read a decimal 1 ulp off on Apple
Silicon (see [r/README.md](r/README.md)). Results never depend on input order.

**k above the number of distinct values.** k is clamped to the number of distinct values (or deck
groups), as Ckmeans.1d.dp does, and `Result.k` reports the bins used. Python also emits a
`UserWarning` when it clamps, as R does; C++ and Scala stay silent and only report k.

## Repository layout

```
cpp/        library (include/lex/lex.hpp, src/), tests/, tools/lex_cli.cpp, cmake/ (fp flags)
python/     nanobind bindings, pyproject (scikit-build-core), pytest suite
scala/      sbt project, munit tests incl. the C++ cross-check, examples (runMain)
r/          R package (Rcpp binding, testthat tests); tools/sync_core.sh copies cpp/ into r/src
examples/   quickstart.py, embeddings.py, twelve.py, cpp/ (CMake consumer)
fixtures/   data/ (committed synthetic data), gen_fixtures.py, make_expected.py, cases.txt, expected_cpp.txt
docs/       algorithm.md
bench/      plot_timing.py
BUILDING.md build and test steps per OS
```

Fixture data is synthetic and committed in `fixtures/data/`; `python fixtures/gen_fixtures.py` checks
it against `fixtures/SHA256SUMS`. `--update` regenerates it, but only on Linux does that reproduce the
committed files bit for bit, because Python's random normals use the platform math library. `python fixtures/make_expected.py
path/to/lex_cli` refreshes `expected_cpp.txt`. The Python tests check that this file is current.

## License

MIT. See [LICENSE](LICENSE).
