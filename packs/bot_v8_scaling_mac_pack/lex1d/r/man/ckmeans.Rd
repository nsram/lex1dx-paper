\name{ckmeans}
\alias{ckmeans}
\alias{lex_hybrid}
\alias{lex_plain}
\alias{ckmeans_exact}
\alias{print.lex1d}
\title{Least-Squares Partition of One-Dimensional Data into k Clusters}
\description{
Partitions \code{x} into \code{k} clusters of contiguous values with a low (by default
near-minimal) within-cluster sum of squares: one-dimensional k-means. The interface and
result follow \code{Ckmeans.1d.dp::Ckmeans.1d.dp} for a fixed \code{k}.

\code{lex_hybrid} (the default method) seeds LEX refinement with an exact dynamic program
over at most \code{coarse_groups} groups. \code{lex_plain} is plain LEX: a greedy Collapse
start, then refinement. \code{ckmeans_exact} is the global optimum by dynamic programming, as
\code{Ckmeans.1d.dp}. \code{ckmeans} is \code{lex_hybrid} with its defaults, a drop-in for
\code{Ckmeans.1d.dp::Ckmeans.1d.dp} with a fixed \code{k}.
}
\usage{
ckmeans(x, k, y = NULL)
lex_hybrid(x, k, y = NULL, coarse_groups = 4096, precision = 0, multi_start = FALSE,
           tolerance = 1e-11)
lex_plain(x, k, y = NULL, precision = 0, multi_start = FALSE, tolerance = 1e-11)
ckmeans_exact(x, k, y = NULL)

\method{print}{lex1d}(x, ...)
}
\arguments{
\item{x}{a numeric vector without \code{NA}, \code{NaN} or infinite values. For
  \code{print}, a \code{"lex1d"} object.}
\item{k}{the number of clusters, a single whole number \eqn{\ge 1}. A \code{k} above the
  number of distinct values (or of precision groups) is reduced to that number with a warning.}
\item{y}{optional non-negative finite weights, one per value of \code{x}, with a positive sum.
  \code{NULL} or a single number means unweighted, as in \code{Ckmeans.1d.dp}.}
\item{coarse_groups}{the number of groups the dynamic program runs on in \code{lex_hybrid}.
  A value \eqn{\le} \code{k} (including 0) gives the \code{lex_plain} result.}
\item{precision}{0 groups exact duplicates only. A positive value also groups values within
  \code{precision} of the running group mean, which changes the problem being solved.}
\item{multi_start}{also try equal-frequency and equal-width starts and refine the best.
  The 2024 deck preset is \code{lex_plain(x, k, precision = 1e-4, multi_start = TRUE)}.}
\item{tolerance}{relative tolerance for a refinement move.}
\item{...}{passed to \code{print} for the numeric components.}
}
\details{
LEX sorts the values, groups equal values (never split across clusters), seeds the clusters,
then moves single groups across cluster boundaries, three clusters at a time, while that lowers
the sum of squares. The \code{lex_hybrid} seed is an exact dynamic program over at most 4096
groups, so for inputs with at most 4096 distinct values its result equals the exact one.
At this default the max SSE gap to exact was under 0.1\% across uniform, normal and
lognormal at \eqn{N = 10^5} and \eqn{10^6} for \eqn{k} from 8 to 512, and under 0.001\%
for \eqn{k \le 64}. The algorithm
is described in \file{docs/algorithm.md} in the source repository.

Among partitions with equal sum of squares every method returns the one with the largest cuts,
as \code{Ckmeans.1d.dp} does.

Differences from \code{Ckmeans.1d.dp}: \code{k} must be a single number (ranges with BIC
selection are not supported and give an error); there are no \code{method} or \code{estimate.k}
arguments and no \code{BIC} component; a \code{y} whose length is neither 1 nor
\code{length(x)} is an error; equal values always share a cluster, whereas
\code{Ckmeans.1d.dp} may split copies of a value that have weight 0 (the sum of squares is the
same).

The core is C++ compiled without floating-point contraction, and the results are bit-identical
to the C++ library and the Python package of the same version for the same input doubles.
Numbers read from text with R's own parser (\code{scan}, \code{read.csv}, \code{as.numeric})
can differ from C \code{strtod} in the last bit on platforms where \code{long double} is
\code{double} (Apple Silicon), so results from the same text file can then differ slightly from
other platforms and from the Python package.
}
\value{
An object of class \code{"lex1d"}, a list with the components of a \code{Ckmeans.1d.dp} result
(except \code{BIC}):
\item{cluster}{integer cluster index (1 to \code{k_used}) of each value of \code{x}, in input order.}
\item{centers}{cluster means (weighted means with \code{y}).}
\item{withinss}{(weighted) within-cluster sum of squares of each cluster.}
\item{size}{number of values per cluster, or the sum of their weights with \code{y} (as in
  \code{Ckmeans.1d.dp}).}
\item{totss}{(weighted) total sum of squares of \code{x}.}
\item{tot.withinss}{\code{sum(withinss)}.}
\item{betweenss}{\code{totss - tot.withinss}.}
\item{xname, yname}{the deparsed \code{x} and \code{y} arguments.}
and these:
\item{count}{integer number of values per cluster.}
\item{upper}{largest value in each cluster.}
\item{sse}{the (weighted) sum of squares as computed by the C++ core; equal to
  \code{tot.withinss} up to rounding.}
\item{k, k_used}{the requested \code{k} and the number of clusters used.}
\item{method}{\code{"lex_hybrid"}, \code{"lex_plain"} or \code{"ckmeans_exact"}.}
\item{start}{which start refinement began from: \code{"coarse_exact"}, \code{"collapse"},
  \code{"equal_frequency"}, \code{"equal_width"} or \code{"exact"}.}
\item{atoms}{number of distinct values (or precision groups).}
\item{sweeps, moves}{refinement sweeps run and moves applied.}
\item{presorted}{\code{TRUE} if \code{x} was already sorted, so the sort was skipped.}
}
\references{
Fisher, W. D. (1958). On grouping for maximum homogeneity. \emph{Journal of the American
Statistical Association}, 53(284), 789--798. \doi{10.1080/01621459.1958.10501479}

Wang, H. and Song, M. (2011). Ckmeans.1d.dp: optimal k-means clustering in one dimension by
dynamic programming. \emph{The R Journal}, 3(2), 29--33. \doi{10.32614/RJ-2011-015}
}
\examples{
set.seed(1)
x <- c(rnorm(500), rnorm(300, mean = 4))
r <- ckmeans(x, 2)
r$centers
r$size
table(r$cluster)

# exact optimum, plain LEX, the deck preset, weights
ckmeans_exact(x, 8)$tot.withinss
lex_plain(x, 8)$tot.withinss
lex_plain(x, 8, precision = 1e-4, multi_start = TRUE)$tot.withinss
lex_hybrid(x, 8, y = rexp(length(x)))$size
}
