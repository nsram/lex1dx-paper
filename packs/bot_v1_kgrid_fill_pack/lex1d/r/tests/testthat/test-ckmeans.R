test_that("every method matches a brute-force optimum or stays above it on small inputs", {
  set.seed(1)
  for (trial in 1:60) {
    n <- sample(5:30, 1)
    x <- round(rnorm(n), sample(c(1, 3), 1)) # rounding gives ties
    k <- sample(1:6, 1)
    best <- brute_sse(x, k)
    e <- suppressWarnings(ckmeans_exact(x, k))
    d <- suppressWarnings(lex_hybrid(x, k))
    expect_equal(e$sse, best, tolerance = 1e-10)
    expect_equal(d$sse, best, tolerance = 1e-10) # n <= 4096: the coarse start is the full DP
    expect_gte(suppressWarnings(lex_plain(x, k))$sse, best - 1e-12)
    expect_gte(suppressWarnings(lex_plain(x, k, precision = 1e-4, multi_start = TRUE))$sse, best - 1e-12)
    expect_identical(suppressWarnings(ckmeans(x, k)), d) # ckmeans() is lex_hybrid()
  }
})

test_that("result is shaped like Ckmeans.1d.dp's and is internally consistent", {
  set.seed(2)
  x <- rlnorm(500)
  r <- ckmeans(x, 6)
  expect_s3_class(r, "lex1d")
  expect_true(all(c("cluster", "centers", "withinss", "size", "totss", "tot.withinss", "betweenss",
                    "xname", "yname") %in% names(r)))
  expect_type(r$cluster, "integer")
  expect_equal(length(r$cluster), length(x))
  expect_equal(as.vector(table(factor(r$cluster, 1:6))), r$count)
  expect_equal(r$size, as.double(r$count))
  expect_equal(r$centers, as.vector(tapply(x, r$cluster, mean)), tolerance = 1e-14)
  expect_equal(r$withinss, as.vector(tapply(x, r$cluster, function(v) sum((v - mean(v))^2))), tolerance = 1e-12)
  expect_equal(r$tot.withinss, r$sse, tolerance = 1e-12)
  expect_equal(r$totss, sum((x - mean(x))^2))
  expect_equal(r$betweenss, r$totss - r$tot.withinss)
  expect_equal(r$k_used, 6L)
  expect_false(r$presorted)
  expect_true(ckmeans(sort(x), 6)$presorted)
  expect_identical(ckmeans(rev(x), 6)$cluster, rev(r$cluster)) # input order does not matter
  expect_output(print(r), "returns 6 clusters")
})

test_that("ties go to the largest cuts, as in Ckmeans.1d.dp", {
  for (f in list(ckmeans, lex_hybrid, lex_plain, ckmeans_exact)) expect_equal(f(c(0, 1, 2), 2)$count, c(2L, 1L))
  expect_equal(lex_plain(0:9, 3)$count, c(4L, 3L, 3L))
  expect_equal(ckmeans_exact(0:9, 3)$count, c(4L, 3L, 3L))
})

test_that("with atoms <= coarse_groups, lex_hybrid is the exact optimum and Refine makes 0 moves", {
  set.seed(4096)
  g <- rnorm(4096)
  sets <- list(
    ties = as.double(1:8),
    small_d = c(0, 0, 0, 1, 1, 2, 2, 2, 2, 5, 9, 9),
    integers = as.double(sample(0:39, 3000, replace = TRUE)),
    gauss300 = rnorm(300),
    gauss4096 = c(g, g) # d = 4096 = the default coarse_groups
  )
  for (x in sets) {
    d <- length(unique(x))
    for (k in c(2, 3, 7, 16, 64)) {
      if (k >= d) next
      e <- ckmeans_exact(x, k)
      for (m in unique(c(d, d + 1, 4096))) {
        if (m < d) next
        r <- lex_hybrid(x, k, coarse_groups = m)
        expect_identical(r$start, "coarse_exact")
        expect_identical(r$count, e$count)
        expect_identical(r$sse, e$sse)
        expect_identical(r$moves, 0L)
      }
      # atoms = coarse_groups + 1: one Collapse merge, so only SSE >= exact is guaranteed
      expect_gte(lex_hybrid(x, k, coarse_groups = d - 1)$sse, e$sse * (1 - 1e-12))
    }
  }
})

test_that("input too large for double arithmetic is an error that names the limit", {
  limit <- sqrt(.Machine$double.xmax) / 4
  # = 3.351951982485649e153; written as 2^510 * (1 - 2^-53) because R on macOS arm64 parses that
  # decimal literal 1 ulp high (as 2^510)
  expect_identical(limit, 2^510 * (1 - 2^-53))
  lim <- "exceeds the limit sqrt(DBL_MAX) / 4 = 3.351951982485649e+153"
  a <- limit / 2 # {-a, a} has mean 0, so n * max |x - mean| = 2a = limit
  up <- a * (1 + .Machine$double.eps)
  for (f in list(ckmeans, lex_plain, ckmeans_exact)) {
    expect_true(is.finite(f(c(-a, a), 2)$tot.withinss))
    expect_error(f(c(-up, up), 2), paste("n * max |x - mean| =", "3.3519519824856493e+153", lim), fixed = TRUE)
  }
  expect_error(ckmeans(rep(c(-2^505, 2^505), each = 500), 2),
               paste(lim, "(n = 1000, max |x - mean| = 1.0474849945267654e+152)"), fixed = TRUE)
  expect_error(ckmeans(c(1e308, 1e308), 1), "mean is not finite")
  b <- limit / 8 # W = 8
  expect_true(is.finite(ckmeans(c(-b, b), 2, y = c(4, 4))$tot.withinss))
  expect_error(ckmeans(c(-b, b) * (1 + .Machine$double.eps), 2, y = c(4, 4)), "max(W, 1) * max |x - mean| =", fixed = TRUE)
  expect_error(ckmeans(c(0, 1), 2, y = c(limit, limit)), paste("the total weight W = 6.703903964971298e+153", lim), fixed = TRUE)
})

test_that("bad input gives clear errors", {
  expect_error(ckmeans(numeric(0), 1), "empty")
  expect_error(ckmeans(c(1, NA), 1), "NA or NaN")
  expect_error(ckmeans(c(1, NaN), 1), "NA or NaN")
  expect_error(ckmeans(c(1, Inf), 1), "Inf")
  expect_error(ckmeans(c(1, -Inf), 1), "Inf")
  expect_error(ckmeans("a", 1), "numeric vector")
  expect_error(ckmeans(matrix(1:4, 2), 1), "numeric vector")
  expect_error(ckmeans(1:5, 0), "whole number >= 1")
  expect_error(ckmeans(1:5, 2.5), "whole number >= 1")
  expect_error(ckmeans(1:5, NA), "whole number >= 1")
  expect_error(ckmeans(1:5, c(1, 9)), "BIC selection")
  expect_error(ckmeans(1:5, 2, y = 1:3), "same length")
  expect_error(ckmeans(1:3, 2, y = c(1, -1, 1)), ">= 0")
  expect_error(ckmeans(1:3, 2, y = c(1, NA, 1)), "finite")
  expect_error(ckmeans(1:3, 2, y = c(0, 0, 0)), "positive sum")
  expect_error(ckmeans(1:3, 2, preset = "plain"), "unused argument")
  expect_error(ckmeans_exact(1:3, 2, coarse_groups = 10), "unused argument")
  expect_error(lex_plain(1:3, 2, coarse_groups = 10), "unused argument")
  expect_error(lex_hybrid(1:3, 2, coarse_groups = -1), "coarse_groups")
  expect_error(lex_plain(1:3, 2, multi_start = NA), "multi_start")
  expect_error(lex_plain(1:3, 2, precision = -1), "precision")
})

test_that("k = 1, constant input and k above the number of distinct values", {
  x <- c(3, 1, 1, 2, 3, 3, 2)
  r1 <- ckmeans(x, 1)
  expect_equal(r1$cluster, rep(1L, 7))
  expect_equal(r1$sse, sum((x - mean(x))^2), tolerance = 1e-15)
  expect_warning(r <- ckmeans(c(5, 5, 5), 2), "k = 2 exceeds the number of distinct values \\(1\\); using k = 1")
  expect_equal(r$sse, 0)
  expect_equal(r$k_used, 1L)
  expect_output(print(r), "returns 1 clusters")
  w <- c(1, 2, 1, 2, 1, 1, 1)
  for (call in list(list(ckmeans), list(lex_plain), list(ckmeans_exact), list(lex_hybrid, y = w))) {
    f <- call[[1]]
    args <- call[-1]
    expect_warning(r <- do.call(f, c(list(x, 10), args)), "using k = 3")
    expect_no_warning(r3 <- do.call(f, c(list(x, 3), args)))
    expect_equal(r$count, c(2L, 2L, 3L))
    expect_identical(r$cluster, r3$cluster)
    expect_equal(r$k, 10L)
    expect_equal(r$k_used, 3L)
  }
  expect_warning(lex_plain(x, 10, precision = 1e-4, multi_start = TRUE), "groups at precision 0.0001 \\(3\\)")
  expect_warning(lex_hybrid(c(1, 1.05, 2), 3, precision = 0.1), "groups at precision 0.1 \\(2\\)")
})

test_that("weights, including zeros, match a weighted brute force", {
  set.seed(3)
  for (trial in 1:40) {
    n <- sample(6:25, 1)
    x <- round(rnorm(n), 2)
    w <- sample(c(0, 0.5, 1, 3), n, replace = TRUE)
    if (sum(w) == 0) w[1] <- 1
    k <- sample(1:5, 1)
    e <- suppressWarnings(ckmeans_exact(x, k, y = w))
    expect_equal(e$sse, brute_sse(x, k, w), tolerance = 1e-10)
    expect_equal(e$size, as.vector(tapply(w, factor(e$cluster, seq_len(e$k_used)), sum)))
  }
  x <- rnorm(300)
  expect_equal(ckmeans(x, 8, y = rep(1, 300))$sse, ckmeans(x, 8)$sse, tolerance = 1e-12)
  expect_identical(ckmeans(x, 8, y = 1)$sse, ckmeans(x, 8)$sse) # a single y is ignored, as in Ckmeans.1d.dp
})
