package examples

import lex.{Lex, Result}

/** LEX quickstart in Scala: three contrasting distributions side by side. Generates its own data.
  *
  *   cd scala && sbt "runMain examples.Quickstart"
  *
  * Same data generator as examples/cpp/quickstart.cpp (splitmix64 + Box-Muller, N = 200,000,
  * 10 % of values rounded to 2 decimals), k = 16. Times include JIT warm-up on the first rows.
  */
object Quickstart:
  final class Rng(private var s: Long): // splitmix64
    def next(): Long =
      s += 0x9e3779b97f4a7c15L
      var z = s
      z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L
      z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL
      z ^ (z >>> 31)
    def uniform(): Double = (next() >>> 11).toDouble * 1.1102230246251565e-16 // 2^-53
    def normal(): Double =
      val u = 1.0 - uniform()
      val v = uniform()
      math.sqrt(-2.0 * math.log(u)) * math.cos(6.283185307179586 * v)

  def sample(dist: String, rng: Rng, n: Int): Array[Double] =
    Array.fill(n) {
      val v = dist match
        case "uniform" => rng.uniform()
        case "normal"  => rng.normal()
        case _         => math.exp(rng.normal())
      if rng.uniform() < 0.10 then math.rint(v * 100) / 100 else v
    }

  def uniformMinMaxMse(x: Array[Double], k: Int): Double =
    val (lo, hi) = (x.min, x.max)
    val step = (hi - lo) / k
    x.map { v =>
      val q = lo + (math.min(math.floor((v - lo) / step), k - 1.0) + 0.5) * step
      (v - q) * (v - q)
    }.sum / x.length

  def timed(f: => Result): (Result, Double) =
    val t = System.nanoTime()
    val r = f
    (r, (System.nanoTime() - t) / 1e6)

  def main(args: Array[String]): Unit =
    val (n, k) = (200000, 16)
    val rng = Rng(2026)
    println(f"${"dist"}%-10s ${"N"}%8s ${"k"}%3s ${"gap hybrid"}%12s ${"gap plain"}%10s ${"uniform/LEX"}%12s ${"t hybrid"}%10s ${"t plain"}%8s ${"t exact"}%8s")
    for dist <- Seq("uniform", "normal", "lognormal") do
      val x = sample(dist, rng, n)
      val (d, td) = timed(Lex.lexHybrid(x, k))
      val (p, tp) = timed(Lex.lexPlain(x, k))
      val (e, te) = timed(Lex.ckmeansExact(x, k))
      val ratio = uniformMinMaxMse(x, k) / (d.sse / n)
      println(f"$dist%-10s $n%8d ${d.k}%3d ${100 * (d.sse / e.sse - 1)}%11.2e%% ${100 * (p.sse / e.sse - 1)}%9.3f%% $ratio%11.2fx $td%8.1fms $tp%6.1fms $te%6.1fms")
    val x = sample("lognormal", rng, n)
    val w = Array.fill(n)(math.exp(rng.normal()))
    val (rw, ew) = (Lex.lexHybrid(x, w, k), Lex.ckmeansExact(x, w, k))
    println(f"weighted lognormal: gap hybrid ${100 * (rw.sse / ew.sse - 1)}%.2e%%, k used ${rw.k}")
