package lex

/** The tie rule: among equal-SSE partitions, the largest cuts, as Ckmeans.1d.dp (R). */
class TiesTest extends munit.FunSuite:
  val presets = Seq(Options(), Methods.plain, Methods.deck, Options(coarseGroups = 6))

  test("{0,1,2}, k = 2: sizes 2|1 for every preset and ckmeansExact (R gives 2|1)") {
    val xs = Array(0.0, 1.0, 2.0)
    assertEquals(Lex.ckmeansExact(xs, 2).counts.toSeq, Seq(2, 1))
    for o <- presets do assertEquals(Methods.solve(xs, 2, o).counts.toSeq, Seq(2, 1), o.toString)
  }

  test("0..9, k = 3: 4|3|3 (R); 0..4, k = 2..4; duplicates") {
    val ten = Array.tabulate(10)(_.toDouble)
    assertEquals(Lex.ckmeansExact(ten, 3).counts.toSeq, Seq(4, 3, 3))
    for o <- presets do assertEquals(Methods.solve(ten, 3, o).counts.toSeq, Seq(4, 3, 3), o.toString)
    val five = Array.tabulate(5)(_.toDouble)
    assertEquals(Lex.ckmeansExact(five, 2).counts.toSeq, Seq(3, 2))
    assertEquals(Lex.ckmeansExact(five, 3).counts.toSeq, Seq(2, 2, 1))
    assertEquals(Lex.ckmeansExact(five, 4).counts.toSeq, Seq(2, 1, 1, 1))
    assertEquals(Lex.ckmeansExact(Array(2.0, 2, 0, 1, 1, 0), 2).counts.toSeq, Seq(4, 2))
  }

  test("DP keeps the largest cut; Refine prefers larger bounds on equal deltas") {
    val a = Atoms(Array(0.0, 1, 2, 3, 4), 0)
    assertEquals(Exact(a, 4).toSeq, Seq(0, 2, 3, 4, 5))
    val b = Array(0, 1, 4, 5)
    val r = Refine(a, b, 1e-11)
    r.run()
    assertEquals(b.toSeq, Seq(0, 2, 4, 5))
  }

  test("results do not depend on input order") {
    val rnd = scala.util.Random(5)
    for _ <- 0 until 200 do
      val xs = Array.fill(5 + rnd.nextInt(80))(rnd.nextInt(12).toDouble)
      val k = 1 + rnd.nextInt(math.min(6, xs.distinct.length))
      val ys = rnd.shuffle(xs.toSeq).toArray
      for o <- presets do // k above the groups is clamped, never an error
        assertEquals(Methods.solve(ys, k, o).counts.toSeq, Methods.solve(xs, k, o).counts.toSeq)
      assertEquals(Lex.ckmeansExact(ys, k).counts.toSeq, Lex.ckmeansExact(xs, k).counts.toSeq)
  }

  test("Refine's unweighted tie window is needed: grid ties end on the largest cuts") {
    // 0..99 ten times, deck preset (lexPlain), k = 8: several partitions have SSE exactly 13000, but
    // the tie move's computed delta is a few ulps off 0 (a tie only within TieEps * mag). With the
    // window set to 0 the result stays on 130|120|130|... cuts.
    val grid = Array.tabulate(1000)(i => (i % 100).toDouble)
    val want = Seq(130, 130, 130, 130, 120, 120, 120, 120)
    assertEquals(Methods.solve(grid, 8, Methods.deck).counts.toSeq, want)
    assertEquals(Lex.ckmeansExact(grid, 8).counts.toSeq, want)
    assert(Refine.TieEps > 0)
  }
