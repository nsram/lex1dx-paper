package lex

import scala.util.Random

class ExactTest extends munit.FunSuite:
  test("ckmeansExact equals the O(k d^2) DP oracle on random small instances") {
    val rnd = Random(99)
    for _ <- 0 until 3000 do
      val xs = Array.fill(1 + rnd.nextInt(40))(rnd.nextInt(25).toDouble - 12)
      val k = 1 + rnd.nextInt(math.min(8, xs.distinct.length))
      val r = Lex.ckmeansExact(xs, k)
      val (opt, _) = Oracle.optimum(xs, k)
      assertEqualsDouble(r.sse, opt, 1e-9 * (1 + opt), s"${xs.toSeq} k=$k")
      assertEquals(r.start, Start.Exact)
  }

  test("ckmeansExact matches the oracle on rounded and heavy-tailed generated data, k = 2..12") {
    val rnd = Random(2379)
    val times = Array.fill(2400)(math.round((10.2 + 0.15 * rnd.nextGaussian()) * 100) / 100.0)
    val pops = Array.fill(800)(math.round(math.exp(10 + 2 * rnd.nextGaussian())).toDouble)
    for xs <- Seq(times, pops); k <- 2 to 12 do
      val opt = Oracle.optimum(xs, k)._1
      assertEqualsDouble(Lex.ckmeansExact(xs, k).sse, opt, 1e-9 * opt)
  }

  test("ckmeansExact validates like lexHybrid and is shift invariant") {
    assertEquals(Lex.ckmeansExact(Array(1.0, 1.0), 2).k, 1) // k clamped to the atoms, as lexHybrid
    intercept[IllegalArgumentException](Lex.ckmeansExact(Array(Double.NaN), 1))
    val rnd = Random(3)
    val base = Array.fill(5000)(rnd.nextGaussian())
    val ref = Lex.ckmeansExact(base, 16)
    for shift <- Seq(1e3, 1e6, 1e8) do
      val r = Lex.ckmeansExact(base.map(_ + shift), 16)
      assert(math.abs(r.sse - ref.sse) / ref.sse < 1e-6, s"shift $shift: ${r.sse} vs ${ref.sse}")
    assert(Methods.solve(base, 16, Methods.plain).sse >= ref.sse * (1 - 1e-12))
  }
