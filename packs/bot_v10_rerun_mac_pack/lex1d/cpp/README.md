# LEX C++ core

C++17 library: `include/lex/lex.hpp` (public API), `src/` (atoms, Collapse, Refine, exact DP).
No dependencies. The code is compiled without FMA contraction (`-ffp-contract=off` on GCC/Clang,
`/fp:precise` on MSVC 2022; see `cmake/LexFloatingPoint.cmake`), which keeps results bit-identical
across x86-64 and arm64 and matches the Python bindings and the Scala reference.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

- `LEX_BUILD_TESTS` and `LEX_BUILD_CLI` default to ON when this is the top-level project and OFF
  under `add_subdirectory`. Link with `lex::lex`.
- **Tests.** `test_edge`, `test_exact`, `test_quality`, `test_coarse`, `test_weighted` and
  `test_ties` use only generated data. They print a few `ORACLE|…` report lines (how often plain
  LEX hits the optimum on small random instances).
- **Optional deck check.** One test reproduces the 2024 deck's printed MSEs. Set `LEX_SIM` to the
  deck's simulation input directory to run it; that data is not shipped, and the test is skipped
  otherwise.
- **`lex_cli`.**
  `lex_cli FILE K [--preset hybrid|plain|deck|exact] [--coarse M] [--weights WFILE]` fits one
  text file: one `value` or `value,weight` per line, with `#` comments, blank lines, an optional
  header line and CRLF endings handled (the exact rules are in `tools/lex_cli.cpp` and the
  repository README). It prints cuts, sizes, SSE, the hex bit patterns of SSE and centres, and
  counters. The presets map to `lex_hybrid` (the default; `--coarse` sets its group count),
  `lex_plain`, `lex_plain` with `deck()`, and `ckmeans_exact`. `fixtures/make_expected.py` uses
  it, and so do the Python and Scala cross-checks.

The API and the algorithm are documented in `include/lex/lex.hpp` and `../docs/algorithm.md`.
