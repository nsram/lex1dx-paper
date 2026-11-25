# lex1d for R

R package for LEX, a fast least-squares partition of one-dimensional data into k contiguous
clusters (1-D k-means). It calls the same C++ core as the Python package through Rcpp, and its
results are bit-identical to the C++ library. Three methods:

- `lex_hybrid()`: the default. An exact dynamic program on at most 4,096 coarse groups, then LEX
  refinement.
- `lex_plain()`: plain LEX, a greedy Collapse start, then refinement.
- `ckmeans_exact()`: the global optimum by dynamic programming, as Ckmeans.1d.dp.

`ckmeans(x, k, y)` is `lex_hybrid()` with its defaults. It takes the arguments and returns the
components of `Ckmeans.1d.dp::Ckmeans.1d.dp()` for a fixed k, so it can stand in for it.

Not on CRAN yet.

## Install from a checkout

Needs R (tested with 4.5.3), Rcpp, and a C++17 compiler: on Windows Rtools, on macOS the
Xcode command line tools.

```sh
git clone https://github.com/nsram/lex1d && cd lex1d
Rscript -e 'install.packages("Rcpp")'
R CMD INSTALL r
```

`R CMD INSTALL r` runs `r/configure`, which copies the C++ core from `cpp/` into `r/src/` (those
copies are not tracked by git). To build a source tarball, copy the core first, because
`R CMD build` does not run `configure`:

```sh
cd r && sh tools/sync_core.sh && cd ..
R CMD build r
R CMD check --as-cran lex1d_0.1.0.tar.gz
```

Tests (testthat; `Ckmeans.1d.dp` is optional and only used for the agreement test):

```sh
Rscript -e 'install.packages(c("testthat", "Ckmeans.1d.dp"))'
cd r/tests/testthat && Rscript -e 'testthat::test_dir(".", package = "lex1d", load_package = "installed")'
```

The bit-identity and `Ckmeans.1d.dp` agreement tests read the repository's `fixtures/` directory.
They find it from `r/tests/testthat` in a checkout; anywhere else (for example under
`R CMD check`), set `LEX1D_FIXTURES` to that directory, or they are skipped.

## Quick start

```r
library(lex1d)
set.seed(2026)
n <- 200000
samples <- list(uniform = runif(n), normal = rnorm(n), lognormal = rlnorm(n))
for (name in names(samples)) {
  x <- samples[[name]]
  tie <- runif(n) < 0.10
  x[tie] <- round(x[tie], 2)                       # some exact duplicates
  d <- lex_hybrid(x, 16)                           # the default, same as ckmeans(x, 16)
  p <- lex_plain(x, 16)
  e <- ckmeans_exact(x, 16)
  cat(sprintf("%-10s gap hybrid %.2e %%  gap plain %.3f %%  start %s\n", name,
              100 * (d$sse / e$sse - 1), 100 * (p$sse / e$sse - 1), d$start))
}
w <- rlnorm(n)
r <- lex_hybrid(samples$normal, 16, y = w)         # weighted
round(r$centers, 3)
```

Output (R 4.5.3, Linux x86-64):

```
uniform    gap hybrid 2.77e-06 %  gap plain 0.477 %  start coarse_exact
normal     gap hybrid 0.00e+00 %  gap plain 0.036 %  start coarse_exact
lognormal  gap hybrid 3.68e-07 %  gap plain 1.041 %  start coarse_exact
 [1] -2.834 -2.153 -1.706 -1.335 -1.014 -0.722 -0.446 -0.182  0.079  0.340
[11]  0.609  0.894  1.208  1.560  2.001  2.640
```

## Interface

```r
lex_hybrid(x, k, y = NULL, coarse_groups = 4096, precision = 0, multi_start = FALSE,
           tolerance = 1e-11)
lex_plain(x, k, y = NULL, precision = 0, multi_start = FALSE, tolerance = 1e-11)
ckmeans_exact(x, k, y = NULL)
ckmeans(x, k, y = NULL)                  # lex_hybrid(x, k, y)
```

| Argument | Meaning |
|---|---|
| `x` | numeric vector; `NA`, `NaN`, `Inf` and empty input are errors |
| `k` | a single whole number ≥ 1; a k above the number of distinct values is reduced with a warning, as in Ckmeans.1d.dp |
| `y` | optional weights (≥ 0, finite, positive sum, same length as `x`); `NULL` or a single number means unweighted |
| `coarse_groups` | `lex_hybrid` only: groups the dynamic program runs on; a value ≤ k (including 0) gives the `lex_plain` result |
| `precision` | 0 groups exact duplicates only; a positive value also groups values within `precision` of the running group mean |
| `multi_start` | also try equal-frequency and equal-width starts |
| `tolerance` | relative tolerance for a refinement move |

The 2024 deck preset (1e-4 grid, three starts) is `lex_plain(x, k, precision = 1e-4, multi_start = TRUE)`.
The arguments match the Python `lex.lex_hybrid`, `lex.lex_plain` and `lex.ckmeans_exact`.

The result has class `"lex1d"` and these components:

- From Ckmeans.1d.dp: `cluster` (1-based, in input order), `centers`, `withinss`, `size` (sum of
  weights with `y`), `totss`, `tot.withinss`, `betweenss`, `xname`, `yname`.
- Added by lex1d: `count` (values per cluster), `upper` (largest value per cluster), `sse` (from
  the C++ core), `k`, `k_used`, `method`, `start`, `atoms`, `sweeps`, `moves`, `presorted`.

Differences from `Ckmeans.1d.dp()`:

- `k` must be a single number. Ranges such as `k = c(1, 9)` with BIC selection give an error, and
  there is no `BIC` component.
- There are no `method` or `estimate.k` arguments.
- A `y` whose length is neither 1 nor `length(x)` is an error.
- Equal values always share a cluster. Ckmeans.1d.dp can split copies of a value that have weight
  0 across two clusters; the sum of squares is the same.
- `ckmeans()` runs `lex_hybrid()`, not the exact solver. Use `ckmeans_exact()` for the global
  optimum. For inputs with at most 4,096 distinct values the `lex_hybrid()` result is the exact one.
  Above that, on continuous simulated data with up to 1 million points, the relative excess SSE
  stayed below 0.002 % for k ≤ 64 and below 0.08 % for k ≤ 512 (see the main README).

On the 88 fixture cases (`fixtures/cases.txt`, 16 of them weighted) `ckmeans_exact()` matched
Ckmeans.1d.dp 4.3.6 in total within-cluster SS (to 1e-9 relative), `size`, and `cluster` for every
value with positive weight.

## Reading numbers on Apple Silicon

R reads decimal text with its own parser (`scan()`, `read.csv()`, `as.numeric()`), which builds
the number in C `long double`. On Apple Silicon `long double` is the same as `double`, so some
values come out 1 ulp (one step in the last bit) away from the correctly rounded value that C
`strtod`, Python and the C++ tools produce. The fit is still bit-identical for the same doubles,
but the same text file can give results that differ in the last bits from other platforms and from
the Python package. The package's own tests read the fixtures with C `strtod` for this reason.
