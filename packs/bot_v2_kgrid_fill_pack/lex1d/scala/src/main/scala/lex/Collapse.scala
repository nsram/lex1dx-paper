package lex

// Mirrors cpp/src/collapse.cpp.

/** Greedy adjacent merging (Ward / PNN restricted to neighbours): repeatedly merge the adjacent
  * pair whose merge raises SSE least, until `k` bins remain. Equal costs merge the leftmost pair:
  * the heap orders by (cost, id), and a tie means exactly equal doubles (no tolerance). This rule
  * only fixes the start. Refine runs afterwards, and its tie sweeps move bounds to the largest cuts
  * (as R's Ckmeans.1d.dp does), so a leftmost choice here does not decide where final ties land.
  * O(d log d) with an indexed 4-ary heap keyed by each bin's cost of merging into its left neighbour.
  * Returns element bounds: bin i is `b(i) until b(i + 1)`.
  */
private[lex] object Collapse:
  def apply(a: Atoms, k: Int): Array[Int] =
    val d = a.d
    val prev = Array.tabulate(d)(_ - 1)
    val next = Array.tabulate(d)(_ + 1) // d marks the end
    def cost(x: Int): Double = a.mergeCost(a.start(prev(x)), a.start(x), a.start(next(x)))

    val heap = MinHeap(d)
    var id = 1
    while id < d do { heap.append(id, cost(id)); id += 1 }
    heap.heapify()
    var bins = d
    while bins > k do
      bins -= 1
      val x = heap.pop() // merge bin x into its left neighbour p
      val p = prev(x)
      val q = next(x)
      next(p) = q
      if q < d then prev(q) = p
      if p > 0 then heap.update(p, cost(p))
      if q < d then heap.update(q, cost(q))

    val b = new Array[Int](k + 1)
    var bin = 0
    var j = 0
    while bin < d do
      b(j) = a.start(bin)
      j += 1
      bin = next(bin)
    b(k) = a.n
    b
end Collapse

/** 4-ary min-heap over ids 0 until capacity, ordered by (key, id), with key updates (as the C++
  * MinHeap). Keys and ids are stored by heap position, so a node's 4 children are adjacent in
  * memory. The pop order is the (key, id) minimum, so it does not depend on the arity.
  */
private final class MinHeap(capacity: Int):
  private val key = new Array[Double](capacity) // key of the id at heap position i
  private val heap = new Array[Int](capacity)   // id at heap position i
  private val pos = Array.fill(capacity)(-1)     // heap position of id
  private var size = 0

  /** Adds an entry without restoring heap order; call heapify() after the last one. */
  def append(id: Int, k: Double): Unit =
    put(id, k, size)
    size += 1

  def heapify(): Unit =
    if size >= 2 then // nothing to order otherwise
      var i = (size - 2) / 4
      while i >= 0 do { down(i); i -= 1 }

  def pop(): Int =
    val top = heap(0)
    size -= 1
    if size > 0 then
      put(heap(size), key(size), 0)
      down(0)
    pos(top) = -1
    top

  def update(id: Int, k: Double): Unit =
    val i = pos(id)
    val oldK = key(i)
    key(i) = k
    if k < oldK then up(i)
    else if k > oldK then down(i)
    else if k != k || oldK != oldK then { up(i); down(pos(id)) } // NaN key: both directions

  private def put(id: Int, k: Double, i: Int): Unit = { heap(i) = id; key(i) = k; pos(id) = i }
  private def less(ka: Double, a: Int, kb: Double, b: Int) = ka < kb || (ka == kb && a < b)

  private def up(i0: Int): Unit =
    var i = i0
    val id = heap(i)
    val k = key(i)
    var done = false
    while !done && i > 0 do
      val p = (i - 1) >> 2
      if less(k, id, key(p), heap(p)) then
        put(heap(p), key(p), i)
        i = p
      else done = true
    put(id, k, i)

  private def down(i0: Int): Unit =
    var i = i0
    val id = heap(i)
    val k = key(i)
    var done = false
    while !done do
      val l = (i << 2) + 1
      if l >= size then done = true
      else
        var c = l
        var ck = key(l)
        var cid = heap(l)
        val end = if l + 4 < size then l + 4 else size
        var j = l + 1
        while j < end do
          if less(key(j), heap(j), ck, cid) then { c = j; ck = key(j); cid = heap(j) }
          j += 1
        if less(ck, cid, k, id) then
          put(cid, ck, i)
          i = c
        else done = true
    put(id, k, i)
end MinHeap
