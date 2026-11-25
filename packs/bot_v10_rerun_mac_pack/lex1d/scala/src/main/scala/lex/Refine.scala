package lex

// Mirrors cpp/src/refine.cpp.

/** 3-bin-window refinement. For the window (A, X, B) around bin X, each of X's two bounds can
  * stay, give X's edge piece to the neighbour, or take the neighbour's edge piece: 3 x 3 - 1 = 8 moves.
  * A piece is the part of one atom inside one bin, so moves never split an atom further.
  *
  * A sweep visits windows from right (X = k-2) to left (X = 1) and applies each window's best move
  * immediately if it lowers SSE; among exactly equal deltas the move with the larger new bounds wins
  * (right bound first). Sweeps repeat until none applies. With k = 2 the single window is
  * (A = 0, X = 1). Bounds are updated in place.
  *
  * Ties (as Ckmeans.1d.dp: largest cuts). After convergence, tie sweeps apply moves that shift the
  * window's bounds right and whose |delta| <= TieEps * mag (equal SSE up to rounding; exactly 0 with
  * weights), taking the largest new bounds; then gain sweeps run again, at most TieRounds times. Tie sweeps are not
  * counted in `sweeps`; their moves are counted in `moves`.
  *
  * Termination: a gain move is applied only if its computed delta is below -tolerance * mag, where mag
  * bounds the size of the delta's terms. Rounding error is far below tolerance * mag, so every applied
  * gain move strictly lowers the true SSE; no partition can repeat, and there are finitely many. Each
  * tie sweep run strictly raises the bound vector (compared from the last bound), so it ends too.
  *
  * The rounding argument fails when weights span more than double's ~16 digits (say 1e12 and 1e-12 in
  * one row): range sums then lose the small weights and gains can come out with the wrong sign. As a
  * guard, a pass also stops when a sweep ends on a partition already seen in that pass. Without a
  * cycle the guard never fires, so it changes no result.
  */
private[lex] final class Refine(a: Atoms, b: Array[Int], tolerance: Double, onMove: Array[Int] => Unit = _ => ()):
  var sweeps = 0
  var moves = 0
  private var mag = 0.0 // accumulated by shift()
  private val tieEps = if a.weighted then 0.0 else Refine.TieEps // weighted: only exactly zero deltas are ties
  // per window: the first shift's delta and mag for each left op (see window())
  private val firstDelta = new Array[Double](3)
  private val firstMag = new Array[Double](3)

  def run(): Unit =
    val k = b.length - 1
    if k >= 2 then
      val top = if k == 2 then 1 else k - 2
      /** Sweeps until none applies; true if any move was applied. */
      def pass(tie: Boolean): Boolean =
        var any = false
        var improved = true
        val seen = scala.collection.mutable.HashSet.empty[Vector[Int]] // cycle guard
        while improved do
          if !tie then sweeps += 1
          improved = false
          var x = top
          while x >= 1 do
            if window(x, tie) then { improved = true; any = true }
            x -= 1
          if improved && !seen.add(b.toVector) then improved = false // cycle: only with extreme weight ratios
        any
      pass(false)
      var round = 0
      while round < Refine.TieRounds && pass(true) && pass(false) do round += 1

  /** Gain mode: applies the best move around bin x if it lowers SSE by more than tolerance * mag.
    * Tie mode: among moves that shift the bounds right with |delta| <= tieEps * mag (TieEps
    * unweighted, 0 weighted; `tolerance` is not used here), applies the one with the largest new bounds.
    */
  private def window(x: Int, tie: Boolean): Boolean =
    val hasB = x + 2 < b.length
    val al = b(x - 1)
    val l = b(x)
    val r = b(x + 1)
    val br = if hasB then b(x + 2) else r
    val id = a.atomId
    val xPieces = id(r - 1) - id(l) + 1
    val aPieces = id(l - 1) - id(al) + 1
    val bPieces = if hasB then id(br - 1) - id(r) + 1 else 0
    // new position of each bound per op: stay, give, take
    val giveL = math.min(a.start(id(l) + 1), r)
    val takeL = math.max(a.start(id(l - 1)), al)
    val giveR = math.max(a.start(id(r - 1)), l)
    val takeR = if hasB then math.min(a.start(id(r) + 1), br) else r

    // The first shift of the usual order depends only on the left op, so its delta and its mag
    // contribution are computed once per left op (mag starts at 0 for each move, so 0 + m1 == m1
    // exactly and results are unchanged). The give-left/take-right move shifts in the other order.
    var firstDone = 0 // bit lo set once firstDelta(lo) / firstMag(lo) hold that left op's shift
    var best = -1
    var bestDelta = 0.0
    var bestMag = 0.0
    var bestL, bestR = 0
    var j = 0
    while j < 8 do
      val lo = Refine.LeftOp(j)
      val ro = Refine.RightOp(j)
      val gives = (if lo == 1 then 1 else 0) + (if ro == 1 then 1 else 0)
      val viable =
        (lo != 2 || aPieces >= 2) && (ro == 0 || hasB) && (ro != 2 || bPieces >= 2) &&
          (lo == 2 || ro == 2 || xPieces > gives)
      if viable then
        val nl = if lo == 0 then l else if lo == 1 then giveL else takeL
        val nr = if ro == 0 then r else if ro == 1 then giveR else takeR
        mag = 0.0
        // apply X's "take" first so X is never empty mid-move
        val delta =
          if lo == 1 && ro == 2 then shift(l, r, nr, br) + shift(al, l, nl, nr)
          else
            if (firstDone & (1 << lo)) == 0 then
              firstDelta(lo) = shift(al, l, nl, r)
              firstMag(lo) = mag
              firstDone |= 1 << lo
            mag = firstMag(lo)
            firstDelta(lo) + shift(nl, r, nr, br)
        val larger = best >= 0 && (nr > bestR || (nr == bestR && nl > bestL))
        if tie then
          if (nr > r || (nr == r && nl > l)) && math.abs(delta) <= tieEps * mag && (best < 0 || larger) then
            best = j; bestL = nl; bestR = nr
        else if best < 0 || delta < bestDelta || (delta == bestDelta && larger) then
          best = j; bestDelta = delta; bestMag = mag; bestL = nl; bestR = nr
      j += 1
    if best >= 0 && (tie || bestDelta < -tolerance * bestMag) then
      b(x) = bestL
      b(x + 1) = bestR
      moves += 1
      onMove(b)
      true
    else false

  /** SSE change from moving the bound between bins `lo until mid` and `mid until hi` to `to`
    * (one piece of k elements with value g). Closed form: k n2/(n2+k) (g-m2)^2 - k n1/(n1-k) (g-m1)^2,
    * where 1 is the bin losing the piece and 2 the bin gaining it.
    */
  private def shift(lo: Int, mid: Int, to: Int, hi: Int): Double =
    if to == mid then 0.0
    else if a.weighted then weightedShift(lo, mid, to, hi)
    else
      val k = math.abs(to - mid).toDouble
      val g = a.value(a.atomId(math.min(to, mid)))
      val right = to > mid // piece moves from the right bin to the left bin
      val n1 = if right then hi - mid else mid - lo // bin losing the piece
      val m1 = if right then a.mean(mid, hi) else a.mean(lo, mid)
      val n2 = if right then mid - lo else hi - mid // bin gaining it
      val m2 = if right then a.mean(lo, mid) else a.mean(mid, hi)
      val w1 = k * n1 / (n1 - k)
      val w2 = k * n2 / (n2 + k)
      val s1 = math.abs(g) + math.abs(m1)
      val s2 = math.abs(g) + math.abs(m2)
      mag += w2 * s2 * s2 + w1 * s1 * s1
      w2 * (g - m2) * (g - m2) - w1 * (g - m1) * (g - m1)

  /** Weighted shift: the same closed form with weights Wp (piece), W1 (losing bin), W2 (gaining
    * bin); a zero-weight piece costs 0. If the piece is the whole losing bin's weight, the delta
    * is taken from the bins' SSEs before and after.
    */
  private def weightedShift(lo: Int, mid: Int, to: Int, hi: Int): Double =
    val p0 = math.min(to, mid)
    val p1 = math.max(to, mid)
    val wp = a.wsum(p0, p1)
    if wp == 0 then 0.0
    else
      val right = to > mid
      val w1 = if right then a.wsum(mid, hi) else a.wsum(lo, mid)
      val w2 = if right then a.wsum(lo, mid) else a.wsum(mid, hi)
      val m1 = if right then a.mean(mid, hi) else a.mean(lo, mid)
      val m2 = if right then a.mean(lo, mid) else a.mean(mid, hi)
      if w1 > wp then
        val g = a.value(a.atomId(p0))
        val ww1 = wp * w1 / (w1 - wp)
        val ww2 = wp * w2 / (w2 + wp)
        val s1 = math.abs(g) + math.abs(m1)
        val s2 = math.abs(g) + math.abs(m2)
        mag += ww2 * s2 * s2 + ww1 * s1 * s1
        ww2 * (g - m2) * (g - m2) - ww1 * (g - m1) * (g - m1)
      else
        val before = a.rangeSse(lo, mid) + a.rangeSse(mid, hi)
        val after = a.rangeSse(lo, to) + a.rangeSse(to, hi)
        mag += math.abs(before) + math.abs(after)
        after - before
end Refine

private[lex] object Refine:
  //                 move: 0  1  2  3  4  5  6  7
  private val LeftOp  = Array(1, 1, 0, 2, 0, 2, 1, 2) // 0 stay, 1 give, 2 take
  private val RightOp = Array(1, 0, 1, 0, 2, 1, 2, 2)
  /** Most gain/tie alternations (ties drift to the largest cuts). */
  val TieRounds = 8
  /** A tie move's |delta| is at most TieEps * mag: equal SSE up to rounding (4 ulp of its terms). */
  val TieEps: Double = 4 * Math.ulp(1.0)
