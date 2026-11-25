
# Shared input checks and the call into the C++ core. method is "lex_hybrid", "lex_plain" or
# "ckmeans_exact"; o holds the LEX options (ignored by ckmeans_exact).
.lex1d <- function(x, k, y, method, o, xname, yname) {
  if (!is.numeric(x) || !is.null(dim(x)) && length(dim(x)) > 1L)
    stop("x must be a numeric vector")
  if (length(x) == 0L) stop("x is empty")
  if (anyNA(x)) stop("x contains NA or NaN")
  if (any(is.infinite(x))) stop("x contains Inf or -Inf")
  if (length(k) != 1L) {
    stop("k must be a single number; ranges such as k = c(1, 9) with BIC selection ",
         "(as in Ckmeans.1d.dp) are not supported")
  }
  if (!is.numeric(k) || !is.finite(k) || k < 1 || k != round(k) || k > .Machine$integer.max)
    stop("k must be a whole number >= 1")
  k <- as.integer(k)

  w <- NULL
  if (!is.null(y) && !(length(y) == 1L && length(x) != 1L)) { # a single y is ignored, as in Ckmeans.1d.dp
    if (!is.numeric(y) || length(y) != length(x)) stop("y (weights) must be numeric with the same length as x")
    if (anyNA(y) || any(!is.finite(y))) stop("y (weights) must be finite and not NA")
    if (any(y < 0)) stop("y (weights) must be >= 0")
    if (sum(y) <= 0) stop("y (weights) must have a positive sum")
    w <- as.double(y)
  }

  if (!is.numeric(o$coarse_groups) || length(o$coarse_groups) != 1L || !is.finite(o$coarse_groups) ||
      o$coarse_groups < 0 || o$coarse_groups != round(o$coarse_groups))
    stop("coarse_groups must be a whole number >= 0")
  if (!is.logical(o$multi_start) || length(o$multi_start) != 1L || is.na(o$multi_start))
    stop("multi_start must be TRUE or FALSE")

  r <- .lex_fit(as.double(x), k, w, method, as.double(o$precision), o$multi_start,
                as.double(o$tolerance), as.double(o$coarse_groups))
  k_used <- length(r$count)
  if (k_used < k) {
    what <- if (o$precision > 0) sprintf("groups at precision %g", o$precision) else "distinct values"
    warning(sprintf("k = %d exceeds the number of %s (%d); using k = %d", k, what, r$atoms, k_used),
            call. = FALSE)
  }

  x <- as.double(x)
  totss <- if (is.null(w)) sum((x - sum(x) / length(x))^2) else sum(w * (x - sum(x * w) / sum(w))^2)
  tot.withinss <- sum(r$withinss)
  structure(list(
    cluster = r$cluster, centers = r$centers, withinss = r$withinss,
    size = if (is.null(w)) as.double(r$count) else r$wsum,
    totss = totss, tot.withinss = tot.withinss, betweenss = totss - tot.withinss,
    xname = xname, yname = yname,
    count = r$count, upper = r$upper, sse = r$sse, k = k, k_used = k_used, method = method,
    start = r$start, atoms = r$atoms, sweeps = r$sweeps, moves = r$moves, presorted = r$presorted
  ), class = "lex1d")
}

lex_hybrid <- function(x, k, y = NULL, coarse_groups = 4096, precision = 0, multi_start = FALSE,
                       tolerance = 1e-11) {
  .lex1d(x, k, y, "lex_hybrid",
         list(coarse_groups = coarse_groups, precision = precision, multi_start = multi_start, tolerance = tolerance),
         deparse(substitute(x)), deparse(substitute(y)))
}

lex_plain <- function(x, k, y = NULL, precision = 0, multi_start = FALSE, tolerance = 1e-11) {
  .lex1d(x, k, y, "lex_plain",
         list(coarse_groups = 0, precision = precision, multi_start = multi_start, tolerance = tolerance),
         deparse(substitute(x)), deparse(substitute(y)))
}

ckmeans_exact <- function(x, k, y = NULL) {
  .lex1d(x, k, y, "ckmeans_exact", list(coarse_groups = 0, precision = 0, multi_start = FALSE, tolerance = 1e-11),
         deparse(substitute(x)), deparse(substitute(y)))
}

# Drop-in for Ckmeans.1d.dp::Ckmeans.1d.dp with a fixed k: lex_hybrid with its defaults.
ckmeans <- function(x, k, y = NULL) {
  .lex1d(x, k, y, "lex_hybrid", list(coarse_groups = 4096, precision = 0, multi_start = FALSE, tolerance = 1e-11),
         deparse(substitute(x)), deparse(substitute(y)))
}

print.lex1d <- function(x, ...) {
  cat("\nCluster centers:\n")
  print(x$centers, ...)
  cat("\nCluster index associated with each element:\n")
  print(x$cluster, ...)
  cat("\nWithin-cluster sum of squares:\n")
  print(x$withinss, ...)
  cat(sprintf("lex1d (%s, start %s) returns %d clusters of sizes %s\n", x$method, x$start,
              x$k_used, paste(format(x$size), collapse = ", ")))
  cat("  total within-cluster sum of squares:", x$tot.withinss, "\n")
  cat("  between-cluster sum of squares:", x$betweenss, "\n")
  cat("  total sum of squares of input vector:", x$totss, "\n")
  if (x$totss > 0) cat("  (between-SS / total-SS):", round(1000 * x$betweenss / x$totss) / 10, "%\n")
  cat("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}
