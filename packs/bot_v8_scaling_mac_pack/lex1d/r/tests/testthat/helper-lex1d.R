# The shared fixtures live in the repository (../fixtures), not in the package. Tests that need them
# use LEX1D_FIXTURES if set, else the checkout's fixtures/ directory, and skip when neither exists.
fixtures_dir <- function() {
  d <- Sys.getenv("LEX1D_FIXTURES", "")
  if (!nzchar(d)) d <- file.path("..", "..", "..", "fixtures")
  if (!file.exists(file.path(d, "expected_cpp.txt"))) testthat::skip("repository fixtures not found (set LEX1D_FIXTURES)")
  d
}

# Fixture numbers go through C strtod (the internal .strtod), as in lex_cli and the Python tests. R's
# own parser (scan(), as.numeric()) accumulates digits in long double, which is double on Apple
# Silicon, so there some values come out 1 ulp off and the results are no longer bit-identical.
read_fixture <- function(dir, name) {
  x <- .strtod(readLines(file.path(dir, "data", paste0(name, ".txt"))))
  stopifnot(!anyNA(x))
  x
}

# lex_cli preset (fixtures/expected_cpp.txt column 4) -> the R call it stands for
run_preset <- function(preset, x, k, y = NULL) {
  switch(preset,
         hybrid = lex_hybrid(x, k, y = y),
         plain = lex_plain(x, k, y = y),
         deck = lex_plain(x, k, y = y, precision = 1e-4, multi_start = TRUE),
         exact = ckmeans_exact(x, k, y = y),
         stop("unknown preset ", preset))
}

# 16 hex digits of an IEEE-754 double (as lex_cli prints sse_bits) -> the double itself
from_bits <- function(hex) {
  bytes <- as.raw(strtoi(substring(hex, seq(1, 15, 2), seq(2, 16, 2)), 16L))
  readBin(bytes, "double", n = 1L, size = 8L, endian = "big")
}

# Optimal (weighted) SSE by an O(k d^2) DP over distinct values; equal values never split.
brute_sse <- function(x, k, w = rep(1, length(x))) {
  v <- sort(unique(x))
  wv <- vapply(v, function(u) sum(w[x == u]), 0)
  d <- length(v)
  k <- min(k, d)
  cost <- function(i, j) { # bin v[i..j], SSE around the weighted mean (two-pass)
    ww <- wv[i:j]; s <- sum(ww)
    if (s == 0) return(0)
    m <- sum(ww * v[i:j]) / s
    sum(ww * (v[i:j] - m)^2)
  }
  best <- vapply(seq_len(d), function(j) cost(1, j), 0)
  if (k > 1) for (b in 2:k) {
    nb <- rep(Inf, d)
    for (j in b:d) for (i in b:j) nb[j] <- min(nb[j], best[i - 1] + cost(i, j))
    best <- nb
  }
  best[d]
}
