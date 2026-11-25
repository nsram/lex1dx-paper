package lex

import scala.util.Random

/** Weighted fitting: optimality against an O(k d^2) weighted oracle, consistency with the
  * unweighted solver, zero weights, input-order independence and the coarse start over groups.
  * Bit-identity with the C++ weighted path is checked in CrossCheckTest.
  */
class WeightedTest extends munit.FunSuite:
  /** Optimal weighted SSE over atoms (equal values merged, weights summed), weighted Welford costs. */
  def oracle(xs: Array[Double], ws: Array[Double], k: Int): Double =
    val atoms = xs.zip(ws).groupMapReduce(_._1)(p => (p._2, p._2 * p._1))((a, b) => (a._1 + b._1, a._2 + b._2))
      .toVector.sortBy(_._1).map { case (v, (w, _)) => (v, w) }
    val d = atoms.length
    val cost = Array.ofDim[Double](d, d + 1)
    for i <- 0 until d do
      var (wsum, mean, m2) = (0.0, 0.0, 0.0)
      for j <- i until d do
        val (v, w) = atoms(j)
        if w > 0 then
          val w2 = wsum + w
          val delta = v - mean
          mean += delta * w / w2
          m2 += delta * delta * wsum * w / w2
          wsum = w2
        cost(i)(j + 1) = m2
    val best = Array.fill(k + 1, d + 1)(Double.PositiveInfinity)
    best(0)(0) = 0.0
    for b <- 1 to k; j <- b to d; i <- (b - 1) until j do
      best(b)(j) = math.min(best(b)(j), best(b - 1)(i) + cost(i)(j))
    best(k)(d)

  /** Weighted SSE of the result's bins over the sorted (value, weight) pairs, two passes per bin. */
  def sseOf(xs: Array[Double], ws: Array[Double], r: Result): Double =
    val p = xs.zip(ws).sortBy(t => (t._1, t._2))
    var (start, total) = (0, 0.0)
    for c <- r.counts do
      val seg = p.slice(start, start + c)
      val w = seg.map(_._2).sum
      if w > 0 then
        val m = seg.map(t => t._1 * t._2).sum / w
        total += seg.map(t => t._2 * (t._1 - m) * (t._1 - m)).sum
      start += c
    total

  test("weighted ckmeansExact matches the O(k d^2) oracle; LEX is never below it") {
    val rng = Random(20261115)
    for _ <- 0 until 300 do
      val n = 3 + rng.nextInt(30)
      val xs = Array.fill(n)(if rng.nextInt(3) == 0 then rng.nextInt(5).toDouble else rng.nextGaussian())
      val ws = Array.fill(n)(0.05 + 5 * rng.nextDouble())
      val d = xs.distinct.length
      val k = 1 + rng.nextInt(math.min(6, d))
      val o = oracle(xs, ws, k)
      val e = Lex.ckmeansExact(xs, ws, k)
      assertEqualsDouble(e.sse, o, 1e-9 * (1 + o))
      assertEqualsDouble(e.sse, sseOf(xs, ws, e), 1e-9 * (1 + o))
      assertEquals(e.start, Start.Exact)
      for opt <- Seq(Options(), Methods.plain, Options(multiStart = true, coarseGroups = 0)) do
        val r = Methods.solve(xs, ws, k, opt)
        assert(r.sse >= o - 1e-9 * (1 + o), s"$opt: ${r.sse} < $o")
        assertEqualsDouble(r.sse, sseOf(xs, ws, r), 1e-9 * (1 + o))
      assertEqualsDouble(Methods.solve(xs, ws, k).sse, o, 1e-9 * (1 + o)) // coarse start = the DP here
  }

  test("extreme weight ratios terminate (Refine cycle guard)") {
    // Same input as the C++ regression test: without the guard, Refine cycles between two partitions.
    val xs = Array(1.4153503275453148, 0.45275334575251319, 1.6316353555349987, -1.3244473157460113, -0.80571464886831112, 1.4002463070789695, 0.35778628864082129, -0.18037482442377539, -0.24516893024270886, -0.032828437692017064, 0.15315063327204761, 0.37458003340011492, -0.74123225585593044, -0.97659650125445463)
    val ws = Array(1e-12, 1e-12, 1e12, 1e12, 1e-12, 0.0, 1e-12, 0.0, 1e12, 0.0, 0.0, 0.0, 1e-12, 0.0)
    assertEquals(Lex.ckmeansExact(xs, ws, 4).counts.length, 4)
    assertEquals(Methods.solve(xs, ws, 4).counts.length, 4)
  }

  test("unit weights give the unweighted SSE; weights are summed over equal values") {
    val rng = Random(3)
    val xs = Array.fill(3000)(math.rint(rng.nextGaussian() * 40) / 8)
    val ones = Array.fill(xs.length)(1.0)
    for k <- Seq(2, 8, 32) do
      assertEqualsDouble(Lex.ckmeansExact(xs, ones, k).sse, Lex.ckmeansExact(xs, k).sse, 1e-12 * Lex.ckmeansExact(xs, k).sse)
      assertEqualsDouble(Methods.solve(xs, ones, k).sse, Methods.solve(xs, k).sse, 1e-12 * Methods.solve(xs, k).sse)
    // doubling every value's weight by duplicating it is the same weighted problem
    val small = xs.take(200)
    val twice = Array.fill(small.length)(2.0)
    assertEqualsDouble(Lex.ckmeansExact(small, twice, 6).sse, Lex.ckmeansExact(small ++ small, 6).sse, 1e-10)
  }

  test("zero weights: the exact SSE equals the problem with those points dropped") {
    val rng = Random(21)
    for _ <- 0 until 80 do
      val npos = 6 + rng.nextInt(12)
      val nzero = 1 + rng.nextInt(8)
      val xs = Array.fill(npos + nzero)(rng.nextGaussian())
      val ws = Array.tabulate(npos + nzero)(i => if i < npos then 0.2 + 2.8 * rng.nextDouble() else 0.0)
      val k = 1 + rng.nextInt(math.min(5, npos))
      val full = Lex.ckmeansExact(xs, ws, k)
      val dropped = Lex.ckmeansExact(xs.take(npos), ws.take(npos), k)
      assertEqualsDouble(full.sse, dropped.sse, 1e-9 * (1 + dropped.sse))
      assert(Methods.solve(xs, ws, k).sse >= dropped.sse - 1e-9 * (1 + dropped.sse))
  }

  test("results do not depend on input order (equal values ordered by weight)") {
    val rng = Random(5)
    val xs = Array.fill(400)(rng.nextInt(40).toDouble) ++ Array.fill(200)(rng.nextGaussian())
    val ws = Array.fill(xs.length)(math.exp(rng.nextGaussian()))
    for opt <- Seq(Options(), Methods.plain, Methods.deck) do
      val base = Oracle.fields(Methods.solve(xs, ws, 12, opt))
      for _ <- 0 until 5 do
        val perm = rng.shuffle(xs.indices.toVector).toArray
        assertEquals(Oracle.fields(Methods.solve(perm.map(xs), perm.map(ws), 12, opt)), base)
    val e = Oracle.fields(Lex.ckmeansExact(xs, ws, 12))
    val rev = xs.indices.reverse.toArray
    assertEquals(Oracle.fields(Lex.ckmeansExact(rev.map(xs), rev.map(ws), 12)), e)
  }

  test("coarse start over groups (M < atoms) uses the weighted DP and stays within 1 % of exact") {
    val rng = Random(9)
    val xs = Array.fill(5000)(rng.nextGaussian())
    val ws = Array.fill(xs.length)(math.exp(rng.nextGaussian()))
    for k <- Seq(4, 16) do
      val e = Lex.ckmeansExact(xs, ws, k).sse
      for m <- Seq(64, 512) do
        val r = Methods.solve(xs, ws, k, Options(coarseGroups = m))
        assertEquals(r.start, Start.CoarseExact)
        assert(r.sse >= e * (1 - 1e-12) && r.sse <= e * (1 + 1e-2), s"k=$k m=$m: ${r.sse} vs $e")
  }

  test("invalid weights are rejected with a clear message") {
    def rejects(ws: Array[Double], fragment: String): Unit =
      val e = intercept[IllegalArgumentException](Methods.solve(Array(1.0, 2.0, 3.0), ws, 2))
      assert(e.getMessage.contains(fragment), e.getMessage)
    rejects(Array(1.0, 1.0), "weights size must match input")
    rejects(Array(1.0, -1.0, 1.0), "weight at index 1 is -1.0")
    rejects(Array(1.0, Double.NaN, 1.0), "weight at index 1 is NaN")
    rejects(Array(0.0, 0.0, 0.0), "positive finite total")
    rejects(Array(1.0, Double.MaxValue, Double.MaxValue), "positive finite total")
  }
