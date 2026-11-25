# Changelog

Repository: https://github.com/nsram/lex1d (PyPI name `lex1d`, Python import name `lex`; R package `lex1d`).

## Unreleased

- **Faster Collapse, weighted fits and Refine; results unchanged.** Collapse's merge heap is now
  4-ary, with each node's key and id stored together. On the M3 Max the grouping pass is 9-17%
  faster and lex_hybrid is about 7% faster end to end on the 8.2M-sample audio. The weighted DP
  reads prefix sums compacted to the group bounds (about 2x faster on the coarse DP), the weighted
  `Atoms` no longer builds the unused unweighted prefixes (16 MB less per million points), and
  Refine reuses the first shift of each window (about 25% faster Refine). The weighted-DP and
  Refine figures are from box microbenchmarks; the 16 MB is computed, not measured.
  Bit-identical: all fixture runs and the Scala cross-check are unchanged.
- **Scala mirrors these optimizations** (4-ary heap, compact weighted prefixes and
  back-pointers, Refine's first-shift reuse), so the Scala code follows the C++ algorithm step
  for step. Bit-identical to the old Scala and to C++ on the 356 fixture runs and a 3000-case fuzz.
- **Input too large for double arithmetic is now an error.** Such input used to give NaN, inf, or
  a finite but wrong partition. Each entry point now checks one rule on the input before any solver
  runs: n · max |x − mean| ≤ √DBL_MAX / 4 = 3.35e153 (with weights, max(W, 1) · max |x − mean| and
  the total weight W, each against the same limit). On failure it raises an error giving the limit,
  the observed value, n or W, and max |x − mean|. The solvers are unchanged, so input that passes gets
  the same results as before. Input just above the limit is rejected even where it would have fitted.
  See docs/algorithm.md, "Overflow".
- **Tests.** A fixture case (`fma_60`, k = 4) and a C++ test whose SSE bits change if the core is
  built with FMA contraction (the existing cases missed an x86-64 `-march=haswell -O3` build); a
  C++ and Scala test that fails if Refine's unweighted tie window is 0; overflow tests in all four
  languages.
- **Tests: exact when atoms ≤ M.** In C++, Python, R and Scala, `lex_hybrid` with at most
  `coarse_groups` distinct values returns `ckmeans_exact`'s cuts and SSE with 0 Refine moves (ties,
  small d, d = M, and d = 4096 at the default M); with one atom more than M only SSE ≥ exact is
  checked.
- **README and docs/algorithm.md: current measurements.** Accuracy, speed against Ckmeans.1d.dp,
  memory, the 8.2M-sample audio result and the effect of a larger M, measured on a MacBook Pro
  (Apple M3 Max, 128 GB). Older figures, including those from M = 16384, are removed.
- Comments: Collapse's leftmost tie rule only fixes the start; zero-weight bins take the global
  weighted mean; `-ffp-contract=off` in r/src/Makevars must stay despite the R CMD check warning.

## 0.1.0 (unreleased)

First public release. It packages the C++ core, Python bindings and the Scala reference implementation.

- **Three named methods in every language**, replacing `fit(..., preset=)` / `fit_exact`:
  `lex_hybrid` (the default: exact DP on 4096 coarse groups, then Refine), `lex_plain` (Collapse,
  then Refine) and `ckmeans_exact` (the DP alone). Scala spells them `Lex.lexHybrid`,
  `Lex.lexPlain`, `Lex.ckmeansExact`. The deck preset is `lex_plain` with `deck()` (C++),
  `Options.deck` (Scala) or `precision = 1e-4, multi_start = TRUE` (Python, R). The R `ckmeans()`
  stays as the Ckmeans.1d.dp drop-in and runs `lex_hybrid()`. No aliases for the old names, since
  nothing was released. Results are unchanged (all 352 fixture runs bit-identical); the fixture
  setting `default` is now called `hybrid`.
- **R tests read fixture numbers with C `strtod`.** R's own parser builds numbers in `long double`,
  which is `double` on Apple Silicon, so some decimals came out 1 ulp off and the bit-identity test
  failed on macOS CI. The R README notes this for users reading text data.

- **Fix: Refine could loop forever with extreme weight ratios** (1e12 and 1e-12 in one input). A
  Refine pass now stops when a partition repeats; normal results are unchanged (all fixture runs
  byte-identical). Found by the weighted stress test on Windows CI, reproduced on Linux. The test's
  brute-force reference now computes bin SSE around the mean, since the raw-sum form lost digits there.

- **Coarse exact start on by default** (`lex_hybrid`, `coarse_groups = 4096`). Use `lex_plain` for
  the previous Collapse start. At M = 4096 the max SSE gap to exact was under 0.1% across
  uniform, normal and lognormal at N = 1e5 and 1e6 for k from 8 to 512, and under 0.001% for
  k ≤ 64 (v8). See `docs/algorithm.md`.
- **Tie rule aligned with Ckmeans.1d.dp.** Among equal-SSE partitions, every solver returns the one
  with the largest cuts. `ckmeans_exact` runs Refine's tie pass after the DP.
- **Weighted fitting** (C++, Python, R and Scala): optional weights on all three methods.
- **k above the number of distinct values (or deck groups) is clamped** to that number, as
  Ckmeans.1d.dp does, instead of raising an error. `Result.k` / `k()` is the number of bins used.
  Python emits a `UserWarning` when it clamps, as R does; C++ and Scala stay silent.
- **R package** `lex1d` (`r/`, Rcpp on the C++ core): `lex_hybrid()`, `lex_plain()`, `ckmeans_exact()`, and `ckmeans(x, k, y)`, which takes the arguments and
  returns the components of `Ckmeans.1d.dp()` for a fixed k (no k ranges or BIC). Results are
  bit-identical to the C++ CLI on all 352 fixture runs.
- **Python bindings** (`lex1d`, import `lex`), built with nanobind and scikit-build-core. Results are
  bit-identical to the C++ CLI.
- **Scala reference** ported to the same tie rule, coarse start and weighted path. It is bit-identical
  to C++ on the whole fixture set, weighted cases included.
- **Platforms:** Linux, macOS and Windows. The MSVC build uses `/fp:precise`, which does not contract
  from VS 2022 on (`/fp:strict` on older MSVC). CI runs on all three OSes; the workflow has not run on
  GitHub yet. Build steps are in `BUILDING.md`.
- **lex_cli input:** one `value` or one `value,weight` per line. `#` comments, blank lines, a header
  line and CRLF endings are handled; other trailing text on a line is now an error, where it used
  to be silently ignored.
- **Fixtures** are committed in `fixtures/data/` and checked by `fixtures/gen_fixtures.py` (pass
  `--update` to regenerate on Linux; regeneration is not bit-portable across OS math libraries); `make_expected.sh` became `make_expected.py`. Both work without bash.
- **Examples** (`examples/`): a Python quickstart on uniform, normal and lognormal data; per-row 4-bit
  quantization of a synthetic matrix; a C++ CMake consumer; Scala `runMain` examples.
- **Removed:**
  - the experimental `Options::window` wider refinement windows (no meaningful accuracy gain);
  - benchmark harnesses and experiments (kept outside the repository);
  - tests that needed external datasets.
