package lex

import scala.util.Random

class QualityTest extends munit.FunSuite:
  override val munitTimeout = scala.concurrent.duration.Duration(600, "s")

  test("random small instances vs exact DP optimum (report; plain LEX is a local optimum)") {
    val rnd = Random(20261112)
    for multi <- Seq(false, true) do
      var (cases, optimal, worst) = (0, 0, 0.0)
      for _ <- 0 until 3000 do
        val xs = Array.fill(3 + rnd.nextInt(28))(rnd.nextInt(21).toDouble)
        val k = 1 + rnd.nextInt(math.min(6, xs.distinct.length))
        val r = Methods.solve(xs, k, Methods.plain.copy(multiStart = multi))
        val (opt, _) = Oracle.optimum(xs, k)
        val sse = Oracle.sse(xs, r.counts.toSeq)
        assertEqualsDouble(r.sse, sse, 1e-9 * (1 + sse))
        assert(sse >= opt - 1e-9 * (1 + opt), s"below optimum?! ${xs.toSeq} k=$k")
        cases += 1
        if sse <= opt + 1e-9 * (1 + opt) then optimal += 1 else worst = math.max(worst, (sse - opt) / opt)
        // lexHybrid (coarse start) is the DP on these short inputs
        assertEqualsDouble(Methods.solve(xs, k).sse, opt, 1e-9 * (1 + opt))
      println(f"ORACLE multiStart=$multi: $optimal/$cases optimal (${100.0 * optimal / cases}%.1f%%), worst relative gap ${100 * worst}%.1f%%")
  }

  test("every applied move lowers the exact SSE (tie moves keep it); atoms are never split") {
    val rnd = Random(7)
    def check(xs: Array[Double], k: Int, multi: Boolean): Unit =
      val a = Atoms(xs.sorted, 0)
      val b = if multi then Lex.equalFrequency(a.n, k) else Collapse(a, k)
      var last = a.sse(b)
      var prev = b.clone()
      Refine(a, b, 1e-11, bb => {
        val s = a.sse(bb)
        val right = prev.reverse.toSeq.zip(bb.reverse.toSeq).find((p, q) => p != q).exists((p, q) => q > p)
        assert(s < last || (right && s <= last + 1e-12 * (1 + math.abs(last))), s"$s vs $last")
        last = s; prev = bb.clone()
      }).run()
      if !multi then assert(b.forall(e => e == a.n || a.start.contains(e)), "split an atom")
    for _ <- 0 until 500 do
      val xs = Array.fill(5 + rnd.nextInt(60))(rnd.nextInt(40).toDouble)
      val d = xs.distinct.length
      check(xs, 1 + rnd.nextInt(math.min(8, d)), multi = false)
      check(xs, 1 + rnd.nextInt(math.min(8, d)), multi = true)
    check(Array.fill(20000)(rnd.nextGaussian()), 16, multi = false)
  }

  test("translation invariance: shifts up to 1e8 give the same SSE, no crash or hang") {
    val rnd = Random(42)
    val base = Array.fill(50000)(rnd.nextGaussian())
    for o <- Seq(Options(), Methods.plain) do
      val ref = Methods.solve(base, 16, o)
      for shift <- Seq(1e3, 1e6, 1e8) do
        val r = Methods.solve(base.map(_ + shift), 16, o)
        assert(math.abs(r.sse - ref.sse) / ref.sse < 1e-6, s"shift $shift: ${r.sse} vs ${ref.sse}")
    val grouped = Options(precision = 1e-4) // the 2024 code hung on +1e6 with precision 1e-4
    val (g0, g6) = (Methods.solve(base, 16, grouped), Methods.solve(base.map(_ + 1e6), 16, grouped))
    assert(math.abs(g6.sse - g0.sse) / g0.sse < 1e-6)
  }

  test("coarse start: lexHybrid is the DP for <= 4096 atoms and close to it above") {
    val rnd = Random(16384)
    val small = Array.fill(3000)(rnd.nextGaussian())
    for k <- Seq(2, 8, 32) do
      val (d, e) = (Methods.solve(small, k), Lex.ckmeansExact(small, k))
      assertEquals(d.start, Start.CoarseExact)
      assertEquals(d.counts.toSeq, e.counts.toSeq)
    val big = Array.fill(40000)(rnd.nextGaussian())
    for k <- Seq(8, 64) do
      val (d, e, p) = (Methods.solve(big, k), Lex.ckmeansExact(big, k), Methods.solve(big, k, Methods.plain))
      assert(d.sse >= e.sse * (1 - 1e-12) && d.sse <= e.sse * (1 + 1e-4), s"k=$k ${d.sse} vs ${e.sse}")
      assert(p.sse >= e.sse * (1 - 1e-12))
  }

  test("atoms <= M: lexHybrid is ckmeansExact (same counts and SSE) with 0 refine moves; atoms = M + 1 is not exact") {
    val rnd = Random(4096)
    val g = Array.fill(4096)(rnd.nextGaussian())
    val sets = Seq(
      "evenly spaced ties" -> Array.tabulate(8)(i => (i + 1).toDouble),
      "duplicates, d = 5" -> Array[Double](0, 0, 0, 1, 1, 2, 2, 2, 2, 5, 9, 9),
      "integers 0..39" -> Array.fill(3000)(rnd.nextInt(40).toDouble),
      "gaussian, d = 300" -> Array.fill(300)(rnd.nextGaussian()),
      "gaussian, d = 4096 (default M)" -> (g ++ g)
    )
    for (name, xs) <- sets do
      val d = xs.distinct.length
      for k <- Seq(2, 3, 7, 16, 64) if k < d do
        val e = Lex.ckmeansExact(xs, k)
        for m <- Seq(d, d + 1, Options().coarseGroups).distinct if m >= d do
          val r = Lex.lexHybrid(xs, k, Options(coarseGroups = m))
          val at = s"$name k=$k M=$m"
          assertEquals(r.start, Start.CoarseExact, at)
          assertEquals(r.counts.toSeq, e.counts.toSeq, at)
          assertEquals(r.sse, e.sse, at)
          assertEquals(r.moves, 0, at)
        val c = Lex.lexHybrid(xs, k, Options(coarseGroups = d - 1))
        assert(c.sse >= e.sse * (1 - 1e-12), s"$name k=$k")
  }

  test("2024 deck config reproduces the deck's printed MSE (LEX_SIM; skipped if unset)") {
    val sim = Oracle.simDir
    assume(sim.isDefined, "set LEX_SIM to lex/simulation/input from LEX.zip to run this test")
    for (dist, printed) <- Seq("uniform_continuous" -> "0.00002032", "gaussian_continuous" -> "0.0006447",
        "exponential_continuous" -> "0.0005397") do
      val all = Oracle.column(sim.get.resolve(s"$dist.txt"))
      val r = Methods.solve(all.drop(all.length - 1048576), 64, Methods.deck)
      val digits = printed.length - printed.indexWhere(c => c != '0' && c != '.')
      assertEquals(BigDecimal(r.mse).round(java.math.MathContext(digits)), BigDecimal(printed))
  }
