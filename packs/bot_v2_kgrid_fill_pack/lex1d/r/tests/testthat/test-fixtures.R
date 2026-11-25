test_that("bit-identical to the C++ lex_cli on every fixture run (fixtures/expected_cpp.txt)", {
  dir <- fixtures_dir()
  lines <- readLines(file.path(dir, "expected_cpp.txt"))
  lines <- lines[!startsWith(lines, "#")]
  n <- 0L
  for (line in lines) {
    f <- strsplit(line, "|", fixed = TRUE)[[1]]
    kv <- strsplit(f[-(1:4)], "=", fixed = TRUE)
    exp <- setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
    x <- read_fixture(dir, f[1])
    w <- if (f[3] == "-") NULL else read_fixture(dir, f[3])
    r <- suppressWarnings(run_preset(f[4], x, as.integer(f[2]), y = w))
    got <- c(bounds = paste(c(0L, cumsum(r$count)), collapse = ","), sizes = paste(r$count, collapse = ","),
             start = r$start, atoms = as.character(r$atoms), sweeps = as.character(r$sweeps),
             moves = as.character(r$moves))
    expect_identical(got, exp[names(got)], label = line)
    expect_identical(r$sse, from_bits(exp[["sse_bits"]]), label = line)
    expect_identical(r$centers, vapply(strsplit(exp[["centers_bits"]], ",")[[1]], from_bits, 0, USE.NAMES = FALSE),
                     label = line)
    n <- n + 1L
  }
  expect_gte(n, 352L) # 88 cases x 4 presets
})

test_that("fixture numbers are read exactly (C strtod), not through R's parser", {
  expect_identical(.strtod(c("0.1", " 1.3220097193920135", "-2e-3", "0x1p-2", "abc", "1 2")),
                   c(0.1, 1.3220097193920135, -2e-3, 0.25, NA, NA))
})

test_that("cluster assignments agree with Ckmeans.1d.dp on the fixtures", {
  # Ckmeans.1d.dp can split equal values across two clusters when some copies have weight 0 (the
  # SSE is the same); lex1d keeps equal values together. So assignments are compared where y > 0.
  skip_if_not_installed("Ckmeans.1d.dp")
  dir <- fixtures_dir()
  cases <- read.table(file.path(dir, "cases.txt"), comment.char = "#", col.names = c("name", "k", "w"),
                      stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cases))) {
    x <- read_fixture(dir, cases$name[i])
    k <- cases$k[i]
    w <- if (cases$w[i] == "-") NULL else read_fixture(dir, cases$w[i])
    ours <- suppressWarnings(ckmeans_exact(x, k, y = w))
    ref <- suppressWarnings(if (is.null(w)) Ckmeans.1d.dp::Ckmeans.1d.dp(x, k)
                            else Ckmeans.1d.dp::Ckmeans.1d.dp(x, k, y = w))
    label <- paste(cases$name[i], k, cases$w[i])
    expect_equal(ours$tot.withinss, ref$tot.withinss, tolerance = 1e-9, label = label)
    keep <- if (is.null(w)) TRUE else w > 0
    expect_identical(ours$cluster[keep], ref$cluster[keep], label = label)
    expect_equal(ours$size, ref$size, tolerance = 1e-12, label = label)
  }
})
