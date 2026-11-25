package lex

// Mirrors cpp/src/exact.cpp.

/** Exact 1-D least-squares partition (optimal k-means in one dimension) by dynamic programming over
  * items (atoms, or groups of atoms for the coarse start), as in Fisher (1958). Each layer is filled
  * by divide and conquer on the monotone optimal split point, which holds because the bin SSE cost
  * satisfies the quadrangle inequality (Wu 1991; Grønlund et al. 2017).
  *
  * Time O(k d log d) after the sort; memory O(k d) ints for back-pointers plus O(d) doubles.
  * Bin costs use double-double prefix sums of centred values and squares: Q - S^2/n.
  * Ties: split points are scanned right to left with a strict `<`, so among exactly equal costs the
  * largest split index wins (as Ckmeans.1d.dp does).
  */
private[lex] object Exact:
  /** Optimal bounds over the atoms of `a` (weighted costs if `a` is weighted). */
  def apply(a: Atoms, k: Int): Array[Int] =
    if a.weighted then weighted(a, a.start, k)
    else apply(a.start, a.value, k)

  /** Item i covers elements start(i) until start(i + 1) and has centred mean value(i). Optimal
    * among partitions that never split an item. Returns element bounds.
    */
  def apply(start: Array[Int], value: Array[Double], k: Int): Array[Int] =
    val d = value.length
    // prefix sums over items of count*value and count*value^2, double-double (hi, lo)
    val sHi, sLo, qHi, qLo = new Array[Double](d + 1)
    var j = 0
    while j < d do
      val c = (start(j + 1) - start(j)).toDouble
      val v = value(j)
      twoSum(sHi, sLo, j, c * v)
      twoSum(qHi, qLo, j, c * v * v)
      j += 1

    /** SSE of items i until j (i < j). */
    def cost(i: Int, j: Int): Double =
      val n = (start(j) - start(i)).toDouble
      val s = (sHi(j) - sHi(i)) + (sLo(j) - sLo(i))
      val q = (qHi(j) - qHi(i)) + (qLo(j) - qLo(i))
      math.max(0.0, q - s * s / n)

    var prev = Array.tabulate(d + 1)(j => if j == 0 then 0.0 else cost(0, j)) // 1 bin
    var cur = new Array[Double](d + 1)
    val from = Array.ofDim[Int](k + 1, d + 1) // from(b)(j): start item of the last of b bins covering 0 until j
    // rows 0 and 1 are never written; dropping them (as the weighted path does) made this loop
    // about 8-10% slower (minimum DP times, box JVM, N = 1e5, k = 512), so the old layout stays here

    /** Fill cur(j) for j in [lo, hi], knowing the optimal split lies in [optLo, optHi]. */
    def fill(b: Int, lo: Int, hi: Int, optLo: Int, optHi: Int): Unit =
      if lo <= hi then
        val mid = (lo + hi) >>> 1
        var best = Double.PositiveInfinity
        val last = math.min(optHi, mid - 1)
        var arg = last
        var i = last
        while i >= optLo do // right to left, strict <: the largest split index wins ties
          val c = prev(i) + cost(i, mid)
          if c < best then { best = c; arg = i }
          i -= 1
        cur(mid) = best
        from(b)(mid) = arg
        fill(b, lo, mid - 1, optLo, arg)
        fill(b, mid + 1, hi, arg, optHi)

    for b <- 2 to k do
      fill(b, b, d, b - 1, d - 1)
      val t = prev; prev = cur; cur = t

    backtrack(start, from.drop(2), d, k)

  /** Weighted DP over the groups with element bounds `start` (groups never split), with bin costs
    * from the weighted prefix sums of `a` (Q - S²/W, clamped at 0). Returns element bounds.
    */
  def weighted(a: Atoms, start: Array[Int], k: Int): Array[Int] =
    val d = start.length - 1
    // the prefix values at the d + 1 group bounds, copied once into contiguous arrays (the same
    // doubles, so the same costs as a.rangeSse(start(i), start(j)))
    val wh, wl, sh, sl, qh, ql = new Array[Double](d + 1)
    var g = 0
    while g <= d do
      val e = start(g)
      wh(g) = a.wHi(e); wl(g) = a.wLo(e); sh(g) = a.wsHi(e); sl(g) = a.wsLo(e); qh(g) = a.wqHi(e); ql(g) = a.wqLo(e)
      g += 1
    def cost(i: Int, j: Int): Double =
      val w = (wh(j) - wh(i)) + (wl(j) - wl(i))
      if w == 0 then 0.0
      else
        val s = (sh(j) - sh(i)) + (sl(j) - sl(i))
        val q = (qh(j) - qh(i)) + (ql(j) - ql(i))
        val v = q - s * s / w
        if v > 0 then v else 0.0
    var prev = Array.tabulate(d + 1)(j => if j == 0 then 0.0 else cost(0, j)) // 1 bin
    var cur = new Array[Double](d + 1)
    val from = Array.ofDim[Int](math.max(k - 1, 0), d + 1)

    def fill(pv: Array[Double], cu: Array[Double], fr: Array[Int], lo: Int, hi: Int, optLo: Int, optHi: Int): Unit =
      if lo <= hi then
        val mid = (lo + hi) >>> 1
        var best = Double.PositiveInfinity
        val last = math.min(optHi, mid - 1)
        var arg = last
        var i = last
        while i >= optLo do // right to left, strict <: the largest split index wins ties
          val c = pv(i) + cost(i, mid)
          if c < best then { best = c; arg = i }
          i -= 1
        cu(mid) = best
        fr(mid) = arg
        fill(pv, cu, fr, lo, mid - 1, optLo, arg)
        fill(pv, cu, fr, mid + 1, hi, arg, optHi)

    for b <- 2 to k do
      fill(prev, cur, from(b - 2), b, d, b - 1, d - 1)
      val t = prev; prev = cur; cur = t
    backtrack(start, from, d, k)

  private def backtrack(start: Array[Int], from: Array[Array[Int]], d: Int, k: Int): Array[Int] =
    val bounds = new Array[Int](k + 1)
    bounds(k) = start(d)
    var end = d
    for b <- k to 2 by -1 do
      end = from(b - 2)(end)
      bounds(b - 1) = start(end)
    bounds

  private def twoSum(hi: Array[Double], lo: Array[Double], j: Int, x: Double): Unit = Atoms.twoSum(hi, lo, j, x)
end Exact
