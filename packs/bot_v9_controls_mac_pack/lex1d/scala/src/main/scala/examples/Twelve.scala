package examples

import lex.Lex

/** 12 values, k = 3: same output format as examples/twelve.py and examples/cpp/twelve.cpp.
  *   cd scala && sbt "runMain examples.Twelve"
  */
object Twelve:
  def main(args: Array[String]): Unit =
    val x = Array(5.5, 1.2, 9.8, 2.0, 1.0, 5.0, 10.4, 1.9, 6.1, 1.2, 2.1, 5.5)
    val r = Lex.lexHybrid(x, 3)
    def bits(v: Double) = f"${java.lang.Double.doubleToRawLongBits(v)}%016x"
    println("sizes   " + r.counts.mkString(" "))
    println("centers " + r.means.map(c => f"$c%.6f").mkString(" "))
    println(f"sse     ${r.sse}%.6f  bits ${bits(r.sse)}")
    println("centers bits " + r.means.map(bits).mkString(" "))
