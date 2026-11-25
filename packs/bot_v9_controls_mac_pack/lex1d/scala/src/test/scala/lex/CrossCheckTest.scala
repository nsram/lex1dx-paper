package lex

/** Scala vs the C++ core on the shared fixture set (../fixtures). expected_cpp.txt is written by
  * fixtures/make_expected.py from the C++ lex_cli (the Python tests check it is current).
  *
  * Same algorithm, same arithmetic: for every case (weighted ones included) and the presets
  * hybrid, plain, deck and exact, the cuts, the SSE and centre bit patterns, the start, the atom
  * count and the sweep/move counts must be identical. Cases with k above the number of distinct
  * values check that both sides clamp k the same way.
  */
class CrossCheckTest extends munit.FunSuite:
  override val munitTimeout = scala.concurrent.duration.Duration(600, "s")

  private def hex(x: Double) = f"${java.lang.Double.doubleToRawLongBits(x)}%016x"
  private def startName(s: Start) = s match
    case Start.Collapse       => "collapse"
    case Start.EqualFrequency => "equal_frequency"
    case Start.EqualWidth     => "equal_width"
    case Start.Exact          => "exact"
    case Start.CoarseExact    => "coarse_exact"

  test("bit-identical to the C++ core on every fixture case, weighted and clamped included") {
    val dir = Oracle.fixturesDir
    val data = dir.resolve("data")
    assume(java.nio.file.Files.isRegularFile(data.resolve("micro_3.txt")),
      s"fixture data missing in $data: run python3 fixtures/gen_fixtures.py")
    val lines = scala.io.Source.fromFile(dir.resolve("expected_cpp.txt").toFile).getLines().filterNot(_.startsWith("#")).toVector
    val cache = scala.collection.mutable.Map.empty[String, Array[Double]]
    var (checked, weighted, clamped) = (0, 0, 0)
    for line <- lines do
      val f = line.split('|')
      val (name, k, w, preset) = (f(0), f(1).toInt, f(2), f(3))
      val cpp = f.drop(4).map(_.split("=", 2)).map(a => a(0) -> a(1)).toMap
      locally {
        val xs = cache.getOrElseUpdate(name, Oracle.column(data.resolve(s"$name.txt")))
        val r =
          if w == "-" then preset match
            case "exact"  => Lex.ckmeansExact(xs, k)
            case "plain"  => Lex.lexPlain(xs, k)
            case "deck"   => Lex.lexPlain(xs, k, Methods.deck)
            case "hybrid" => Lex.lexHybrid(xs, k)
          else
            val ws = cache.getOrElseUpdate(w, Oracle.column(data.resolve(s"$w.txt")))
            weighted += 1
            preset match
              case "exact"  => Lex.ckmeansExact(xs, ws, k)
              case "plain"  => Lex.lexPlain(xs, ws, k)
              case "deck"   => Lex.lexPlain(xs, ws, k, Methods.deck)
              case "hybrid" => Lex.lexHybrid(xs, ws, k)
        if r.k < k then clamped += 1
        val ours = Map(
          "bounds" -> r.counts.scanLeft(0)(_ + _).mkString(","),
          "sizes" -> r.counts.mkString(","),
          "sse" -> cpp("sse"), // readable copy; sse_bits is compared
          "sse_bits" -> hex(r.sse),
          "centers_bits" -> r.means.map(hex).mkString(","),
          "start" -> startName(r.start),
          "atoms" -> r.atoms.toString,
          "sweeps" -> r.sweeps.toString,
          "moves" -> r.moves.toString
        )
        assertEquals(ours, cpp, s"$name k=$k w=$w $preset")
        checked += 1
      }
    println(s"CROSSCHECK|scala vs c++|$checked cases identical|$weighted weighted|$clamped with k clamped")
    assert(checked == lines.size && weighted >= 48 && clamped >= 16)
  }
