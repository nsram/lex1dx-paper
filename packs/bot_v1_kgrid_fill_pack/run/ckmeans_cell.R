#!/usr/bin/env Rscript
# One Ckmeans.1d.dp cell in its own process.
# Usage:
#   Rscript ckmeans_cell.R INPUT.f64 K METHOD WARM REPS SORTED BOUNDS_OUT
# INPUT.f64 is little-endian IEEE-754 float64 (raw or already std::sort-ed).
# SORTED is 1 if INPUT is already in nondecreasing order, 0 if not.
# Times only the Ckmeans.1d.dp() call (1 unrecorded warm-up + REPS timed
# reps). The median uses the same upper-middle index as v10 (0-based length/2).
# Ckmeans.1d.dp always sorts internally; when SORTED=1 the call is on a
# pre-sorted array, so that internal sort is not excluded.
# Writes k+1 0-based bin edges (left-to-right on the number line) to BOUNDS_OUT
# and prints one CELL line on stdout.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: Rscript ckmeans_cell.R INPUT.f64 K METHOD WARM REPS SORTED BOUNDS_OUT")
}
path <- args[[1]]
k <- as.integer(args[[2]])
method <- args[[3]]
warm <- as.integer(args[[4]])
reps <- as.integer(args[[5]])
is_sorted <- as.integer(args[[6]])
bounds_out <- args[[7]]
if (is.na(k) || k < 1L || is.na(warm) || warm < 0L || is.na(reps) || reps < 1L) {
  stop("K, WARM, REPS must be integers; K >= 1, REPS >= 1, WARM >= 0")
}
if (!method %in% c("linear", "loglinear", "quadratic")) {
  stop("METHOD must be linear, loglinear, or quadratic")
}

info <- file.info(path)
if (is.na(info$size) || info$size %% 8 != 0) stop(paste("bad f64 file:", path))
n <- as.integer(info$size / 8)
con <- file(path, "rb")
x <- readBin(con, what = "double", n = n, size = 8L, endian = "little")
close(con)
if (length(x) != n) stop(paste("short read:", path))
if (any(!is.finite(x))) stop(paste("non-finite value in", path))

suppressPackageStartupMessages(library(Ckmeans.1d.dp))
ver <- as.character(packageVersion("Ckmeans.1d.dp"))

# Distinct count on this vector (for the skip_k_ge_D guard). For a pre-sorted
# input this matches v10's D after the +0 / std::sort canonicalisation.
if (is_sorted == 1L) {
  D <- if (n == 0L) 0L else 1L + sum(x[-1] != x[-n])
} else {
  D <- length(unique(x))
}
if (k > D || as.integer(k) >= n) {
  cat(sprintf(
    "CELL method=Ckmeans.1d.dp k=%d n=%d D=%d status=skip_k_ge_D ckmeans_method=%s version=%s sorted=%d\n",
    k, n, D, method, ver, is_sorted
  ))
  quit(save = "no", status = 0)
}

now <- function() proc.time()[["elapsed"]]
res <- NULL
if (warm > 0L) {
  for (i in seq_len(warm)) {
    res <- Ckmeans.1d.dp(x, k = k, method = method)
  }
}
times <- numeric(reps)
for (i in seq_len(reps)) {
  t0 <- now()
  res <- Ckmeans.1d.dp(x, k = k, method = method)
  times[[i]] <- now() - t0
}
st <- sort(times)
# Same median as v10 scale_cell: 0-based index length/2, upper middle for even n.
med <- st[[as.integer(length(st) / 2) + 1L]]

# Bin edges on the number line: clusters ordered by increasing centre.
if (is.null(res) || is.null(res$centers) || is.null(res$size)) {
  stop("Ckmeans.1d.dp returned no centres / sizes")
}
o <- order(res$centers)
sizes <- as.integer(res$size[o])
if (sum(sizes) != n) stop(sprintf("cluster sizes sum to %d, not n=%d", sum(sizes), n))
bounds <- c(0L, cumsum(sizes))
writeLines(as.character(bounds), bounds_out)

tot <- if (!is.null(res$tot.withinss)) res$tot.withinss else NA_real_
k_eff <- length(sizes)
cat(sprintf(
  paste0(
    "CELL method=Ckmeans.1d.dp k=%d n=%d D=%d status=ok ",
    "ckmeans_method=%s version=%s sorted=%d ",
    "post_s=%.9f time_min_s=%.9f time_max_s=%.9f ",
    "tot_withinss=%.17g k_eff=%d warm=%d reps=%d\n"
  ),
  k, n, D, method, ver, is_sorted,
  med, min(times), max(times),
  tot, k_eff, warm, reps
))
