package lex

class EdgeCasesTest extends munit.FunSuite:
  def rejects(xs: Array[Double], k: Int, fragment: String, o: Options = Options()): Unit =
    val e = intercept[IllegalArgumentException](Methods.solve(xs, k, o))
    assert(e.getMessage.contains(fragment), e.getMessage)

  test("invalid input is rejected with a clear message") {
    rejects(Array.empty, 1, "empty")
    rejects(Array(1.0, Double.NaN), 1, "index 1 is NaN")
    rejects(Array(Double.PositiveInfinity), 1, "index 0 is Infinity")
    rejects(Array(0.0, Double.NegativeInfinity), 1, "index 1 is -Infinity")
    rejects(Array(1.0, 2.0), 0, "k must be at least 1")
    rejects(Array(1.0), 1, "precision", Options(precision = -1))
    rejects(Array(1.0), 1, "tolerance", Options(tolerance = 0))
    rejects(Array(1.0), 1, "coarseGroups", Options(coarseGroups = -1))
  }

  test("input too large for double arithmetic is rejected: n * max |x - mean| <= sqrt(DBL_MAX) / 4") {
    val L = math.sqrt(Double.MaxValue) / 4
    assertEquals(L, 3.351951982485649e153)
    val lim = "exceeds the limit sqrt(DBL_MAX) / 4 = 3.351951982485649E153"
    def all(xs: Array[Double], k: Int, fragment: String): Unit =
      for o <- Seq(Options(), Methods.plain, Methods.deck) do rejects(xs, k, fragment, o)
      val e = intercept[IllegalArgumentException](Lex.ckmeansExact(xs, k))
      assert(e.getMessage.contains(fragment), e.getMessage)
    def passes(xs: Array[Double], k: Int): Unit =
      for o <- Seq(Options(), Methods.plain, Methods.deck) do assert(Methods.solve(xs, k, o).sse.isFinite)
      assert(Lex.ckmeansExact(xs, k).sse.isFinite)
    // n = 2 at the limit: {-a, a} has mean exactly 0, so n * max |x - mean| = 2a
    val a = L / 2
    val aUp = Math.nextUp(a)
    passes(Array(-a, a), 2)
    all(Array(-aUp, aUp), 2, s"n * max |x - mean| = 3.3519519824856493E153 $lim")
    all(Array(-aUp, aUp), 2, "(n = 2, max |x - mean| = 1.6759759912428246E153)")
    // 1000 values: half at -2^505, half at +2^505 (exact sums, so the mean is 0 on any platform)
    val runs = Array.tabulate(1000)(i => if i < 500 then -math.pow(2, 505) else math.pow(2, 505))
    all(runs, 2, s"$lim (n = 1000, max |x - mean| = 1.0474849945267654E152)")
    // a mean that is not finite: the sum of the values overflows
    all(Array(1e308, 1e308), 1, "the mean is not finite because the sum of the values overflows (limit 1.7976931348623157E308)")
    all(Array(1.7e308, -1e308, 1.7e308, -1e308), 2, "mean is not finite")
    // the case that broke Refine's rounding bound before (cuts 3|1, SSE 1.04e308) is above the limit
    val b = math.sqrt(Double.MaxValue) * 0.93 / 2
    all(Array(-b, -b, b, b * (1 + 1.0 / 256)), 2, lim)
    // data just below the limit match the same data scaled by 2^-600 (exact scaling)
    val c = L / 1000 / 1.2 // mean c / 6, so n * max |x - mean| = 0.97 L
    val edge = Array.tabulate(1000)(i => if i % 3 == 0 then -c else if i % 3 == 1 then c * 0.5 else c)
    val scaled = edge.map(x => java.lang.Math.scalb(x, -600))
    for k <- Seq(2, 3, 5) do
      for o <- Seq(Options(), Methods.plain) do // not deck: its precision does not scale
        val (r, q) = (Methods.solve(edge, k, o), Methods.solve(scaled, k, o))
        assertEquals(r.counts.toSeq, q.counts.toSeq); assertEquals(r.moves, q.moves); assertEquals(r.sweeps, q.sweeps)
        assert(r.sse.isFinite)
      assertEquals(Lex.ckmeansExact(edge, k).counts.toSeq, Lex.ckmeansExact(scaled, k).counts.toSeq)
    // weighted: max(W, 1) * max |x - mean| and W against the same limit
    def wrejects(xs: Array[Double], ws: Array[Double], fragment: String): Unit =
      for f <- Seq(() => Lex.lexHybrid(xs, ws, 2), () => Lex.lexPlain(xs, ws, 2), () => Lex.ckmeansExact(xs, ws, 2)) do
        val e = intercept[IllegalArgumentException](f())
        assert(e.getMessage.contains(fragment), e.getMessage)
    def wpasses(xs: Array[Double], ws: Array[Double]): Unit =
      assert(Lex.lexHybrid(xs, ws, 2).sse.isFinite && Lex.lexPlain(xs, ws, 2).sse.isFinite && Lex.ckmeansExact(xs, ws, 2).sse.isFinite)
    val a4 = L / 8 // W = 8: 8 * a4 = L
    wpasses(Array(-a4, a4), Array(4.0, 4.0))
    wrejects(Array(-Math.nextUp(a4), Math.nextUp(a4)), Array(4.0, 4.0), s"max(W, 1) * max |x - mean| = 3.3519519824856493E153 $lim (W = 8.0, ")
    wpasses(Array(-L, L), Array(0.25, 0.25)) // W < 1 counts as 1
    wrejects(Array(-Math.nextUp(L), Math.nextUp(L)), Array(0.25, 0.25), "(W = 0.5, max |x - mean| = 3.3519519824856493E153)")
    wrejects(Array(0.0, 1.0), Array(L, L), s"the total weight W = 6.703903964971298E153 $lim (max |x - mean| = 0.5)")
    wrejects(Array(-1e200, 0, 1e200), Array(1e10, 1, 1e10), s"max(W, 1) * max |x - mean| = 2.0000000001E210 $lim")
  }

  test("weighted Atoms carry no unweighted prefix sums") {
    val w = Atoms(Array(1.0, 2.0, 2.0), Array(1.0, 0.5, 2.0), 0)
    assert(w.weighted && w.wHi.length == 4)
  }

  test("k above the number of atoms is clamped to it (as Ckmeans.1d.dp); Result.k reports the bins used") {
    val xs = Array(3.0, 1.0, 1.0, 2.0, 3.0, 3.0, 2.0) // 3 distinct values
    for o <- Seq(Options(), Methods.plain, Methods.deck); k <- Seq(4, 7, 1000) do
      val (r, r3) = (Methods.solve(xs, k, o), Methods.solve(xs, 3, o))
      assertEquals(r.k, 3); assertEquals(r.counts.toSeq, Seq(2, 2, 3)); assertEquals(r.sse, 0.0)
      assertEquals(Oracle.fields(r), Oracle.fields(r3))
    assertEquals(Lex.ckmeansExact(xs, 10).k, 3)
    assertEquals(Lex.ckmeansExact(Array(1.0, 1.0), 2).counts.toSeq, Seq(2))
    val g = Methods.solve(Array(1.0, 1.05, 2.0), 3, Options(precision = 0.1)) // precision groups count as atoms
    assertEquals(g.k, 2); assertEquals(g.atoms, 2); assertEquals(g.counts.toSeq, Seq(2, 1))
    val w = Array(1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0)
    assertEquals(Oracle.fields(Methods.solve(xs, w, 5)), Oracle.fields(Methods.solve(xs, w, 3))); assertEquals(Methods.solve(xs, w, 5).k, 3)
    assertEquals(Oracle.fields(Lex.ckmeansExact(xs, w, 5)), Oracle.fields(Lex.ckmeansExact(xs, w, 3)))
  }

  test("all-equal values: k = 1 has SSE 0") {
    val r = Methods.solve(Array.fill(7)(3.5), 1)
    assertEquals(r.counts.toSeq, Seq(7)); assertEquals(r.sse, 0.0); assertEquals(r.means(0), 3.5)
  }

  test("k = 1 gives the total sum of squares") {
    val xs = Array(1.0, 2.0, 4.0, 8.0)
    assertEqualsDouble(Methods.solve(xs, 1).sse, Oracle.sse(xs, Seq(4)), 1e-12)
  }

  test("k = number of distinct values gives SSE 0 and keeps duplicates together") {
    val xs = Array(3.0, 1.0, 2.0, 2.0, 5.0, 1.0)
    val r = Methods.solve(xs, 4)
    assertEquals(r.counts.toSeq, Seq(2, 2, 1, 1)); assertEquals(r.sse, 0.0)
  }

  test("tiny inputs, k = 2 refinement, and signed zeros") {
    assertEquals(Methods.solve(Array(0.0, 1.0), 2).counts.toSeq, Seq(1, 1))
    assertEquals(Methods.solve(Array(0.0, 0.0, 1.0, 2.0), 2).counts.toSeq, Seq(2, 2)) // optimum
    assertEquals(Methods.solve(Array(-0.0, 0.0, 0.0), 1).atoms, 1)
  }

  test("input array is not modified; assign() maps every value to its bin") {
    val xs = Array(9.0, 1.0, 5.0, 2.0, 7.0, 3.0)
    val copy = xs.clone()
    val r = Methods.solve(xs, 3)
    assertEquals(xs.toSeq, copy.toSeq)
    val labels = xs.sorted.map(r.assign)
    assertEquals(labels.toSeq, r.counts.toSeq.zipWithIndex.flatMap((c, i) => Seq.fill(c)(i)))
  }

  test("Collapse ties merge the leftmost pair first") {
    val xs = Array.tabulate(12)(_.toDouble) // every adjacent merge costs the same
    assertEquals(Collapse(Atoms(xs.clone(), 0), 11).toSeq, Seq(0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
  }
