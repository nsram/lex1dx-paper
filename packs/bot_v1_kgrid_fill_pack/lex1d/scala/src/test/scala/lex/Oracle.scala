package lex

/** The tests name each LEX method by its Options: coarseGroups = 0 runs `Lex.lexPlain`, anything
  * else `Lex.lexHybrid` (the cross-check calls both directly).
  */
object Methods:
  val plain: Options = Options(coarseGroups = 0)
  val deck: Options = Options.deck.copy(coarseGroups = 0) // Lex.lexPlain(xs, k, Options.deck)
  def solve(xs: Array[Double], k: Int, o: Options = Options()): Result =
    if o.coarseGroups == 0 then Lex.lexPlain(xs, k, o) else Lex.lexHybrid(xs, k, o)
  def solve(xs: Array[Double], ws: Array[Double], k: Int, o: Options): Result =
    if o.coarseGroups == 0 then Lex.lexPlain(xs, ws, k, o) else Lex.lexHybrid(xs, ws, k, o)
  def solve(xs: Array[Double], ws: Array[Double], k: Int): Result = solve(xs, ws, k, Options())

/** Exact optimum by dynamic programming over atoms (runs of equal values), O(k d^2), for tests.
  * Bin SSEs use weighted Welford updates, so no large sums of squares are differenced.
  */
object Oracle:
  /** (optimal SSE, counts) for k contiguous bins of the sorted values of xs. */
  def optimum(xs: Array[Double], k: Int): (Double, Vector[Int]) =
    val s = xs.sorted
    val vals = s.distinct
    val cnt = vals.map(v => s.count(_ == v))
    val d = vals.length
    val cost = Array.ofDim[Double](d, d + 1) // cost(i)(j): SSE of atoms i until j
    for i <- 0 until d do
      var (n, mean, m2) = (0.0, 0.0, 0.0)
      for j <- i until d do
        val c = cnt(j).toDouble
        val delta = vals(j) - mean
        val n2 = n + c
        mean += delta * c / n2
        m2 += delta * delta * n * c / n2
        n = n2
        cost(i)(j + 1) = m2
    val best = Array.fill(k + 1, d + 1)(Double.PositiveInfinity)
    val from = Array.ofDim[Int](k + 1, d + 1)
    best(0)(0) = 0.0
    for b <- 1 to k; j <- b to d; i <- (b - 1) until j do
      val c = best(b - 1)(i) + cost(i)(j)
      if c < best(b)(j) then { best(b)(j) = c; from(b)(j) = i }
    var (j, counts) = (d, List.empty[Int])
    for b <- k to 1 by -1 do
      val i = from(b)(j)
      counts = cnt.slice(i, j).sum :: counts
      j = i
    (best(k)(d), counts.toVector)

  /** Exact SSE of given counts over sorted xs (two-pass per bin). */
  def sse(xs: Array[Double], counts: Seq[Int]): Double =
    val s = xs.sorted
    var (start, total) = (0, 0.0)
    for c <- counts do
      val seg = s.slice(start, start + c)
      val m = seg.sum / c
      total += seg.map(v => (v - m) * (v - m)).sum
      start += c
    total

  /** Every field of a result as comparable values (Result holds IArrays, which compare by reference);
    * doubles as raw bits, so equality means bit-identical.
    */
  def fields(r: Result): Seq[Any] =
    def bits(x: Double) = java.lang.Double.doubleToRawLongBits(x)
    Seq(r.counts.toSeq, r.means.toSeq.map(bits), r.upper.toSeq.map(bits), bits(r.sse), r.start, r.atoms, r.sweeps, r.moves)

  /** The shared fixture folder: LEX_FIXTURES, else ../fixtures from scala/ (or ./fixtures). */
  def fixturesDir: java.nio.file.Path =
    val ok = (p: java.nio.file.Path) => java.nio.file.Files.isRegularFile(p.resolve("cases.txt"))
    sys.env.get("LEX_FIXTURES").map(java.nio.file.Path.of(_)).filter(ok).getOrElse {
      val cwd = java.nio.file.Path.of(sys.props("user.dir")).toAbsolutePath.normalize
      Seq(cwd.resolve("../fixtures").normalize, cwd.resolve("fixtures")).find(ok).getOrElse(cwd.resolve("../fixtures").normalize)
    }

  def column(path: java.nio.file.Path): Array[Double] =
    val src = scala.io.Source.fromFile(path.toFile)
    try src.getLines().map(_.trim).filter(_.nonEmpty).map(_.toDouble).toArray
    finally src.close()

  /** LEX_SIM = path to LEX.zip's lex/simulation/input (the folder with uniform_continuous.txt). */
  def simDir: Option[java.nio.file.Path] =
    sys.env.get("LEX_SIM").flatMap { raw =>
      val p = java.nio.file.Path.of(raw).toAbsolutePath.normalize
      Seq(p, p.resolve("input"), p.resolve("simulation/input"), p.resolve("lex/simulation/input"))
        .find(c => java.nio.file.Files.isRegularFile(c.resolve("uniform_continuous.txt")))
    }
