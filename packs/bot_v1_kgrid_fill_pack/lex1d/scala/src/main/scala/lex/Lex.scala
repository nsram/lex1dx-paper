package lex

// Mirrors cpp/include/lex/lex.hpp (API) and cpp/src/lex.cpp (validation, starts, summaries).

import java.util.Arrays

/** Which partition refinement started from. */
enum Start:
  case Collapse, EqualFrequency, EqualWidth, Exact, CoarseExact

/** Options for `Lex.lexHybrid` and `Lex.lexPlain` (`Lex.ckmeansExact` takes none). The defaults
  * solve the exact problem (no rounding).
  *
  * @param precision  0 groups exact duplicates only. A positive value groups each sorted value
  *                   with the current run if it lies within `precision` of the run's mean
  *                   (the 2024 deck rule). This changes the problem being solved.
  * @param multiStart also try equal-frequency and equal-width partitions and refine from
  *                   whichever has the lowest SSE (ties prefer collapse, then frequency).
  * @param tolerance  a move is applied only if it lowers SSE by more than `tolerance` times
  *                   the size of the terms in its delta (a bound on rounding error).
  * @param coarseGroups `lexHybrid` only (`lexPlain` ignores it; default 4096). Collapse stops at
  *                   M = max(k, min(coarseGroups, atoms)) groups, the exact DP cuts them into k bins,
  *                   and refinement continues. M = k (coarseGroups <= k, including 0) is
  *                   `lexPlain`'s Collapse start; M = atoms starts from the optimum.
  */
final case class Options(
    precision: Double = 0.0,
    multiStart: Boolean = false,
    tolerance: Double = 1e-11,
    coarseGroups: Int = 4096
)

object Options:
  /** The 2024 deck setting: 1e-4 grouping, three-way start. The deck preset is `Lex.lexPlain(xs, k, Options.deck)`. */
  val deck: Options = Options(precision = 1e-4, multiStart = true)

/** A fitted quantizer. Bin `i` holds the next `counts(i)` sorted values, whose largest value is `upper(i)`. */
final case class Result(
    counts: IArray[Int],
    means: IArray[Double],
    upper: IArray[Double],
    sse: Double,
    start: Start,
    atoms: Int,
    sweeps: Int,
    moves: Int
):
  def k: Int = counts.length
  def n: Int = counts.sum
  def mse: Double = sse / n

  /** Bin index for value `x`: the first bin whose upper value is >= x. */
  def assign(x: Double): Int =
    var (lo, hi) = (0, k - 1)
    while lo < hi do
      val mid = (lo + hi) >>> 1
      if upper(mid) >= x then hi = mid else lo = mid + 1
    lo

object Lex:
  /** The default method: the coarse exact start (Collapse to `options.coarseGroups` groups, the
    * exact DP on them), then LEX refinement. Partitions `xs` into `k` contiguous bins of its sorted
    * values, minimizing the sum of squared deviations from the bin means. The result is a local
    * optimum: no 3-bin-window move improves it. A `k` above the number of atoms (distinct values,
    * or precision groups) is clamped to it, as Ckmeans.1d.dp does; `Result.k` is the number of bins used.
    */
  def lexHybrid(xs: Array[Double], k: Int, options: Options = Options()): Result =
    val (atoms, kk) = prepare(xs, k, options)
    solve(atoms, kk, options)

  /** Weighted `lexHybrid`: `weights` holds a finite weight >= 0 per value (same length, positive
    * total). Minimizes the weighted SSE; `Result.sse` is the weighted SSE and `Result.means` the
    * weighted bin means. Equal values are ordered by weight, so results never depend on input order.
    * A bin whose values all have weight 0 has no weighted mean; its mean is the weighted mean of all
    * the data, and it adds 0 to the SSE.
    */
  def lexHybrid(xs: Array[Double], weights: Array[Double], k: Int, options: Options): Result =
    val (atoms, kk) = prepareWeighted(xs, weights, k, options)
    solve(atoms, kk, options)

  def lexHybrid(xs: Array[Double], weights: Array[Double], k: Int): Result = lexHybrid(xs, weights, k, Options())

  /** Plain LEX: the Collapse start, then refinement (`options.coarseGroups` is ignored). */
  def lexPlain(xs: Array[Double], k: Int, options: Options = Options()): Result =
    lexHybrid(xs, k, options.copy(coarseGroups = 0))

  /** Weighted `lexPlain`. */
  def lexPlain(xs: Array[Double], weights: Array[Double], k: Int, options: Options): Result =
    lexHybrid(xs, weights, k, options.copy(coarseGroups = 0))

  def lexPlain(xs: Array[Double], weights: Array[Double], k: Int): Result = lexPlain(xs, weights, k, Options())

  /** The globally optimal partition, by dynamic programming (see `Exact`), as Ckmeans.1d.dp: O(k d log d) time and
    * O(k d) memory for d distinct values. Use it as a reference or when d * k is moderate
    * (about 1e5 distinct values in well under a second, 1e6 in seconds at k = 64).
    */
  def ckmeansExact(xs: Array[Double], k: Int): Result =
    val (atoms, kk) = prepare(xs, k, Options())
    exactLargestCuts(atoms, kk)

  /** Weighted `ckmeansExact`. */
  def ckmeansExact(xs: Array[Double], weights: Array[Double], k: Int): Result =
    val (atoms, kk) = prepareWeighted(xs, weights, k, Options())
    exactLargestCuts(atoms, kk)

  private def solve(atoms: Atoms, k: Int, options: Options): Result =
    val m = math.max(k, math.min(options.coarseGroups, atoms.d))
    val (first, firstStart) =
      if m <= k then (Collapse(atoms, k), Start.Collapse)
      else if m == atoms.d then (Exact(atoms, k), Start.CoarseExact) // every atom is a group
      else // Collapse to m groups, then the DP over the groups
        val groups = Collapse(atoms, m)
        if atoms.weighted then (Exact.weighted(atoms, groups, k), Start.CoarseExact)
        else
          val means = Array.tabulate(m)(g => atoms.mean(groups(g), groups(g + 1)))
          (Exact(groups, means, k), Start.CoarseExact)
    val (bounds, start) = if options.multiStart then bestStart(atoms, first, firstStart, k) else (first, firstStart)
    val refine = Refine(atoms, bounds, options.tolerance)
    refine.run()
    summarize(atoms, bounds, start, refine.sweeps, refine.moves)

  /** The DP's optimum, then the refine tie pass: ties hidden by the DP's rounding also end on the
    * largest cuts; gain moves cannot apply at the optimum, so only tie moves are counted.
    */
  private def exactLargestCuts(atoms: Atoms, k: Int): Result =
    val bounds = Exact(atoms, k)
    val refine = Refine(atoms, bounds, Options().tolerance)
    refine.run()
    summarize(atoms, bounds, Start.Exact, 0, refine.moves)

  private def checkCommon(xs: Array[Double]): Unit =
    if xs.isEmpty then fail("input is empty")
    val bad = xs.indexWhere(x => !java.lang.Double.isFinite(x))
    if bad >= 0 then fail(s"input value at index $bad is ${xs(bad)}; only finite values are allowed")

  private def checkOptions(k: Int, options: Options): Unit =
    if k < 1 then fail(s"k must be at least 1, got $k")
    if !(options.precision >= 0 && options.precision < Double.PositiveInfinity) then
      fail(s"precision must be finite and >= 0, got ${options.precision}")
    if !(options.tolerance > 0 && options.tolerance < 1) then fail(s"tolerance must be in (0, 1), got ${options.tolerance}")
    if options.coarseGroups < 0 then fail(s"coarseGroups must be >= 0, got ${options.coarseGroups}")

  /** Validates input and returns sorted, grouped data and k clamped to the number of atoms. */
  private def prepare(xs: Array[Double], k: Int, options: Options): (Atoms, Int) =
    checkCommon(xs)
    checkOptions(k, options)
    val sorted = xs.map(_ + 0.0) // -0.0 + 0.0 is +0.0, so -0 and +0 get the same bits
    Arrays.sort(sorted)
    val atoms = Atoms(sorted, options.precision)
    checkRange(atoms)
    (atoms, math.min(k, atoms.d))

  private def prepareWeighted(xs: Array[Double], weights: Array[Double], k: Int, options: Options): (Atoms, Int) =
    checkCommon(xs)
    if weights == null || weights.length != xs.length then fail("weights size must match input")
    var total = 0.0
    var i = 0
    while i < weights.length do
      val w = weights(i)
      if !java.lang.Double.isFinite(w) || w < 0 then fail(s"weight at index $i is $w; weights must be finite and >= 0")
      total += w
      i += 1
    if !(total > 0) || !java.lang.Double.isFinite(total) then fail(s"weights must sum to a positive finite total, got $total")
    checkOptions(k, options)
    // order by (value, weight), so equal values with different weights get one fixed order
    val x = xs.map(_ + 0.0)
    val order = (0 until x.length).toArray.sortWith((p, q) => x(p) < x(q) || (x(p) == x(q) && weights(p) < weights(q)))
    val atoms = Atoms(order.map(x(_)), order.map(weights(_)), options.precision)
    checkRange(atoms)
    (atoms, math.min(k, atoms.d))

  private def fail(msg: String): Nothing = throw IllegalArgumentException(msg)

  /** The overflow rule (docs/algorithm.md, "Overflow"), checked once on the input before any solver
    * runs; the same rule and wording as check_range in cpp/src/lex.cpp, with doubles spelled the
    * Java way. Z = max |x - mean| (the mean the solvers centre on), N = n or the total weight W:
    * max(N, 1) * Z and N must be at most sqrt(Double.MaxValue) / 4. A mean that is not finite
    * gives Z = inf and fails too. Input that passes runs exactly as before.
    */
  private def checkRange(a: Atoms): Unit =
    val limit = math.sqrt(Double.MaxValue) / 4 // 3.351951982485649E153
    val z = math.max(math.abs(a.sorted(0) - a.center), math.abs(a.sorted(a.n - 1) - a.center))
    val total = if a.weighted then a.wsum(0, a.n) else a.n.toDouble
    val product = math.max(total, 1.0) * z
    if !(product <= limit && total <= limit) then
      val lim = s"the limit sqrt(DBL_MAX) / 4 = $limit"
      val zs = s"max |x - mean| = $z"
      val what =
        if !(total <= limit) then s"the total weight W = $total exceeds $lim ($zs)"
        else if a.weighted then s"max(W, 1) * max |x - mean| = $product exceeds $lim (W = $total, $zs)"
        else s"n * max |x - mean| = $product exceeds $lim (n = ${a.n}, $zs)"
      val mean =
        if java.lang.Double.isFinite(a.center) then ""
        else s"; the mean is not finite because the sum of the values overflows (limit ${Double.MaxValue})"
      fail(s"input range too large for double arithmetic: $what$mean")

  private def bestStart(a: Atoms, first: Array[Int], firstStart: Start, k: Int): (Array[Int], Start) =
    val candidates = Seq(
      (first, firstStart),
      (equalFrequency(a.n, k), Start.EqualFrequency),
      (equalWidth(a.sorted, k), Start.EqualWidth)
    ).filter(_._1 != null)
    candidates.minBy(c => a.sse(c._1)) // minBy keeps the first of equal values

  /** First k - n%k bins get n/k values, the rest one more (as in the 2024 code). */
  private[lex] def equalFrequency(n: Int, k: Int): Array[Int] =
    val f = n / k
    val extra = n - f * k
    val b = new Array[Int](k + 1)
    for i <- 1 to k do b(i) = b(i - 1) + f + (if i > k - extra then 1 else 0)
    b

  /** Bins of equal value width; null if any bin would be empty. The last bin always ends at n. */
  private[lex] def equalWidth(sorted: Array[Double], k: Int): Array[Int] =
    val n = sorted.length
    val lo = sorted(0)
    val width = (sorted(n - 1) - lo) / k
    val b = new Array[Int](k + 1)
    for i <- 1 until k do b(i) = upperBound(sorted, lo + i * width)
    b(k) = n
    if (1 to k).forall(i => b(i) > b(i - 1)) then b else null

  /** Number of values <= x. */
  private def upperBound(sorted: Array[Double], x: Double): Int =
    var lo = 0
    var hi = sorted.length
    while lo < hi do
      val mid = (lo + hi) >>> 1
      if sorted(mid) <= x then lo = mid + 1 else hi = mid
    lo

  private def summarize(a: Atoms, b: Array[Int], start: Start, sweeps: Int, moves: Int): Result =
    val k = b.length - 1
    val means = new Array[Double](k)
    var sse = 0.0
    for i <- 0 until k do // per bin: two passes over its values, centred
      val (lo, hi) = (b(i), b(i + 1))
      if a.weighted then
        var w = 0.0
        var s = 0.0
        var e = lo
        while e < hi do { w += a.elemWeight(e); s += a.elemWeight(e) * (a.sorted(e) - a.center); e += 1 }
        // a bin whose values all have weight 0 has no weighted mean: its mean falls back to the
        // global weighted mean (centred 0) and it adds nothing to the SSE
        val m = if w == 0 then 0.0 else s / w
        if w > 0 then
          e = lo
          while e < hi do
            val d = a.sorted(e) - a.center - m
            sse += a.elemWeight(e) * d * d
            e += 1
        means(i) = a.center + m
      else
        var s = 0.0
        var e = lo
        while e < hi do { s += a.sorted(e) - a.center; e += 1 }
        val m = s / (hi - lo)
        e = lo
        while e < hi do
          val d = a.sorted(e) - a.center - m
          sse += d * d
          e += 1
        means(i) = a.center + m
    Result(
      IArray.tabulate(k)(i => b(i + 1) - b(i)),
      IArray.unsafeFromArray(means),
      IArray.tabulate(k)(i => a.sorted(b(i + 1) - 1)),
      sse,
      start,
      a.d,
      sweeps,
      moves
    )
end Lex
