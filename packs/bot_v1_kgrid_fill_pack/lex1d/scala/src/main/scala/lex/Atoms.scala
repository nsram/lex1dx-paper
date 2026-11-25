package lex

// Mirrors cpp/src/atoms.cpp and the Atoms struct in cpp/src/detail.hpp.

/** Sorted data, centred on its mean and grouped into atoms: runs of values treated as one value.
  * Atom `j` covers elements `start(j) until start(j + 1)`, all with centred value `value(j)`;
  * `atomId(e)` is the atom of element `e`. Bin sums come from double-double prefix sums, so a bin's
  * mean depends only on its bounds (not on the history of moves) and is accurate to a few ulps.
  *
  * Weighted data (`weighted`) also carry a weight per sorted element and double-double prefix sums
  * of W = Σw, S = Σw·z and Q = Σw·z² (z the centred atom value); means, merge costs and SSEs then
  * use the weighted formulas. The unweighted formulas are unchanged. Weighted Atoms have no
  * unweighted prefix sums (`prefixHi` / `prefixLo` are null, as the C++ vectors are empty); `mean`
  * reads the weighted prefixes instead.
  */
private[lex] final class Atoms private (
    val sorted: Array[Double],
    val center: Double,
    val atomId: Array[Int],
    val start: Array[Int],
    val value: Array[Double],
    prefixHi: Array[Double], // null if weighted
    prefixLo: Array[Double], // null if weighted
    val elemWeight: Array[Double], // null if unweighted
    private[lex] val wHi: Array[Double], private[lex] val wLo: Array[Double],
    private[lex] val wsHi: Array[Double], private[lex] val wsLo: Array[Double],
    private[lex] val wqHi: Array[Double], private[lex] val wqLo: Array[Double]
):
  val n: Int = sorted.length
  val d: Int = value.length
  val weighted: Boolean = elemWeight != null

  /** Sum of weights of elements `l until r` (weighted only). */
  def wsum(l: Int, r: Int): Double = (wHi(r) - wHi(l)) + (wLo(r) - wLo(l))
  /** Weighted sum of centred values of elements `l until r` (weighted only). */
  def wS(l: Int, r: Int): Double = (wsHi(r) - wsHi(l)) + (wsLo(r) - wsLo(l))
  /** Weighted sum of squared centred values of elements `l until r` (weighted only). */
  def wQ(l: Int, r: Int): Double = (wqHi(r) - wqHi(l)) + (wqLo(r) - wqLo(l))

  /** Weighted SSE of elements `l until r` about their weighted mean (0 if the weight is 0). */
  def rangeSse(l: Int, r: Int): Double =
    val w = wsum(l, r)
    if w == 0 then 0.0
    else
      val s = wS(l, r)
      val q = wQ(l, r)
      val v = q - s * s / w
      if v > 0 then v else 0.0

  /** Mean centred value of elements `l until r` (l < r); the weighted mean if weighted. A range of
    * total weight 0 has mean 0, i.e. the global weighted mean `center` (a bin of zero-weight values).
    */
  def mean(l: Int, r: Int): Double =
    if weighted then
      val w = wsum(l, r)
      if w == 0 then 0.0 else wS(l, r) / w
    else ((prefixHi(r) - prefixHi(l)) + (prefixLo(r) - prefixLo(l))) / (r - l)

  /** SSE increase from merging bins `l until m` and `m until r`: na nb / (na + nb) (mean_a - mean_b)^2,
    * with total weights for na, nb if weighted (0 if either weight is 0).
    */
  def mergeCost(l: Int, m: Int, r: Int): Double =
    if weighted then
      val w1 = wsum(l, m)
      val w2 = wsum(m, r)
      if w1 == 0 || w2 == 0 then 0.0
      else
        val diff = mean(l, m) - mean(m, r)
        w1 * w2 / (w1 + w2) * diff * diff
    else
      val na = (m - l).toDouble
      val nb = (r - m).toDouble
      val diff = mean(l, m) - mean(m, r)
      na * nb / (na + nb) * diff * diff

  /** SSE of the atom values for bounds `b` (two-pass per bin; prefix sums per bin if weighted). */
  def sse(b: Array[Int]): Double =
    var total = 0.0
    var i = 0
    if weighted then
      while i < b.length - 1 do { total += rangeSse(b(i), b(i + 1)); i += 1 }
      return total
    while i < b.length - 1 do
      val m = mean(b(i), b(i + 1))
      var e = b(i)
      while e < b(i + 1) do
        val dv = value(atomId(e)) - m
        total += dv * dv
        e += 1
      i += 1
    total
end Atoms

private[lex] object Atoms:
  /** Groups sorted values: precision 0 joins exact duplicates; precision > 0 joins a value to the
    * current run if it is less than `precision` above the run's mean (the 2024 deck rule). Equal
    * values always share an atom.
    */
  def apply(sorted: Array[Double], precision: Double): Atoms =
    val n = sorted.length
    var s = 0.0
    var i = 0
    while i < n do { s += sorted(i); i += 1 }
    val m = s / n
    var c = 0.0
    i = 0
    while i < n do { c += sorted(i) - m; i += 1 }
    val center = m + c / n

    val atomId = new Array[Int](n)
    val start = new Array[Int](n + 1)
    val value = new Array[Double](n)
    var d = 0
    i = 0
    while i < n do
      val first = i
      var sum = sorted(i) - center
      atomId(i) = d
      i += 1
      while i < n && (sorted(i) == sorted(i - 1) || (precision > 0 && sorted(i) - center - sum / (i - first) < precision)) do
        sum += sorted(i) - center
        atomId(i) = d
        i += 1
      start(d) = first
      value(d) = if precision == 0 then sorted(first) - center else sum / (i - first)
      d += 1
    start(d) = n

    val (prefixHi, prefixLo) = valuePrefixes(n, value, atomId)
    new Atoms(sorted, center, atomId, java.util.Arrays.copyOf(start, d + 1), java.util.Arrays.copyOf(value, d),
      prefixHi, prefixLo, null, null, null, null, null, null, null)

  /** Same with a non-negative weight per sorted value (aligned with `sorted`; equal values ordered
    * by weight). The precision rule uses the weighted run mean (the plain mean while the run's weight is 0).
    */
  def apply(sorted: Array[Double], weights: Array[Double], precision: Double): Atoms =
    val n = sorted.length
    var wt = 0.0
    var st = 0.0
    var i = 0
    while i < n do { wt += weights(i); st += weights(i) * sorted(i); i += 1 }
    val m = st / wt
    var c = 0.0
    i = 0
    while i < n do { c += weights(i) * (sorted(i) - m); i += 1 }
    val center = m + c / wt

    val atomId = new Array[Int](n)
    val start = new Array[Int](n + 1)
    val value = new Array[Double](n)
    var d = 0
    i = 0
    while i < n do
      val first = i
      var wrun = weights(i)
      var srun = weights(i) * (sorted(i) - center)
      var crun = 1
      var sUnw = sorted(i) - center
      atomId(i) = d
      i += 1
      while i < n && (sorted(i) == sorted(i - 1) ||
          (precision > 0 && (if wrun > 0 then sorted(i) - center - srun / wrun
                             else sorted(i) - center - sUnw / crun) < precision)) do
        wrun += weights(i)
        srun += weights(i) * (sorted(i) - center)
        sUnw += sorted(i) - center
        crun += 1
        atomId(i) = d
        i += 1
      start(d) = first
      value(d) = if precision == 0 then sorted(first) - center else if wrun > 0 then srun / wrun else sUnw / crun
      d += 1
    start(d) = n

    val wHi, wLo, wsHi, wsLo, wqHi, wqLo = new Array[Double](n + 1)
    var e = 0
    while e < n do
      val w = weights(e)
      val z = value(atomId(e))
      twoSum(wHi, wLo, e, w)
      twoSum(wsHi, wsLo, e, w * z)
      twoSum(wqHi, wqLo, e, w * z * z)
      e += 1
    new Atoms(sorted, center, atomId, java.util.Arrays.copyOf(start, d + 1), java.util.Arrays.copyOf(value, d),
      null, null, weights, wHi, wLo, wsHi, wsLo, wqHi, wqLo)

  /** Double-double prefix sums of the centred atom value per element. */
  private def valuePrefixes(n: Int, value: Array[Double], atomId: Array[Int]): (Array[Double], Array[Double]) =
    val prefixHi = new Array[Double](n + 1)
    val prefixLo = new Array[Double](n + 1)
    var hi = 0.0
    var lo = 0.0
    var e = 0
    while e < n do
      val x = value(atomId(e))
      val t = hi + x // TwoSum: t + err == hi + x exactly
      val bv = t - hi
      lo += (hi - (t - bv)) + (x - bv)
      hi = t
      prefixHi(e + 1) = hi
      prefixLo(e + 1) = lo
      e += 1
    (prefixHi, prefixLo)

  /** hi(j + 1) + lo(j + 1) = hi(j) + lo(j) + x, with TwoSum's rounding error folded into lo. */
  private[lex] def twoSum(hi: Array[Double], lo: Array[Double], j: Int, x: Double): Unit =
    val s = hi(j) + x
    val bv = s - hi(j)
    hi(j + 1) = s
    lo(j + 1) = lo(j) + (hi(j) - (s - bv)) + (x - bv)
end Atoms
