# LEX: algorithm notes

## Problem

Given values x₁ … xₙ (optionally with weights wᵢ ≥ 0) and k ≥ 1, find k contiguous bins of the
sorted values that minimize the (weighted) sum of squared deviations from the bin means:

    SSE(b) = Σᵢ Σ_{e ∈ bin i} w_e (x_e − m_i)²,   m_i = Σ w_e x_e / Σ w_e

This is 1-D k-means (optimal scalar quantization). Its optimal bins are always contiguous in
sorted order, so a solution is k + 1 cut indices `0 = b₀ < b₁ < … < b_k = n` into the sorted data.

## Atoms

The data are sorted (an O(n) check skips the sort for already-sorted input) and centred on their
mean (with a second correction pass), which makes the arithmetic translation invariant. Equal values form an
**atom**: a run the algorithm never splits (−0.0 counts as +0.0). The number of atoms is d ≤ n.
A k above d is clamped to d, and `Result.k` (C++ `k()`) reports the number of bins used. Under
the deck preset, d counts the precision groups. Ckmeans.1d.dp does the same: R 4.5.3 with
Ckmeans.1d.dp 4.3.6 returns d clusters and warns, for example k = 4 on {1, 1, 2, 2, 3, 3, 3}
gives 3 clusters of sizes 2, 2, 3. The Python bindings warn as well (`UserWarning`). C++ and
Scala clamp silently and report the bins used in `k()` / `Result.k`.

Bin means come from double-double (TwoSum) prefix sums. A mean therefore depends only on the bounds,
is accurate to a few ulps, and does not drift as moves accumulate. No sum of squares is ever
differenced in Collapse or Refine.

## Collapse (greedy start)

Start with one bin per atom. Repeatedly merge the adjacent pair whose merge raises SSE least,
until k bins remain. The merge cost of bins a, b is the closed form

    Δ = n_a n_b / (n_a + n_b) · (m_a − m_b)²

(with weights, total weights replace counts). An indexed 4-ary min-heap keyed by (cost, index)
holds each bin's cost of merging into its left neighbour, so the d − k merges take O(d log d) in
total (C++ and Scala use the same heap). On equal costs the leftmost pair merges first; "equal"
means exactly equal doubles, with no tolerance. This rule only fixes the start: Refine's tie sweeps
then move ties to the largest cuts (see the tie rule below).

## Refine (3-bin window, 8 moves)

For the window (A, X, B) around bin X, each of X's two bounds can **stay**, **give** X's edge piece
to the neighbour, or **take** the neighbour's edge piece. That is 3 × 3 − 1 = **8 moves**. A piece is the part of one atom inside one bin, so moves never
split an atom further. A move's SSE change is the sum of two boundary shifts, each in closed form:

    δ = k n₂/(n₂ + k) (g − m₂)² − k n₁/(n₁ − k) (g − m₁)²

Here a piece of k elements with value g moves from bin 1 to bin 2. X's "take" is applied first, so
X is never empty mid-move.

A sweep visits windows from X = k − 2 down to X = 1 and applies each window's best move at once.
Sweeps repeat until none applies. With k = 2 there is one window.

**Termination.** A gain move is applied only if its computed δ < −tolerance · mag. Here mag is the
sum of the delta's term sizes, and tolerance is 1e-11, far above the rounding error of δ (a small
multiple of 2⁻⁵³ · mag). So every applied gain move strictly lowers the true SSE, and no partition repeats.
The result is a local optimum: no 3-bin-window move improves it.

## Exact DP (`ckmeans_exact`)

Fisher's (1958) dynamic program over atoms:

    D(b, j) = min_{i<j} D(b − 1, i) + SSE(atoms i … j − 1)

Bin costs Q − S²/n come from double-double prefix sums, clamped at 0. Each layer is filled by
divide and conquer on the monotone optimal split point. This is valid because the bin cost
satisfies the quadrangle inequality (Wu 1991; Grønlund et al. 2017). Time is O(k d log d) after the
sort; memory is O(k d) ints for back-pointers. The O(k n) SMAWK variant is not implemented.

## Coarse exact start (`lex_hybrid`, the default; `coarse_groups = 4096`)

1. Collapse stops at M = max(k, min(coarse_groups, d)) groups.
2. The exact DP cuts those M groups (weighted by their counts) into k bins, optimally among
   partitions that never split a group.
3. Refine continues from there at atom level.

M = k is the plain Collapse start (`lex_plain`). M = d (any input with at most 4096 distinct values) starts from
the global optimum, so the result equals `ckmeans_exact`'s (same cuts and SSE) and Refine makes no
moves; the tests check this in all four languages, including d = M. The cost is O(d log d) for Collapse plus
O(k M log M) for the DP.

**Measured basis for the default** (MacBook Pro, Apple M3 Max, 128 GB; C++ library; the source
data are outside this repository, see the README):

- **Excess error.** On continuous simulated data with up to 1 million points, the relative excess
  SSE over the exact optimum grew roughly as k²: below 0.002 % for k ≤ 64 and below 0.08 % for
  k ≤ 512, with no trend in N.
- **Time.** Faster than Ckmeans.1d.dp in every continuous cell from k = 8, at least 2.2× faster for
  k ≥ 16, and a median 68× faster at N = 1 million, k = 512. On 8.2 million audio samples at
  k = 512: 3.0 s and 0.61 GB of peak memory with 0.028 % excess SSE, against 172 s and 68 GB for
  Ckmeans.1d.dp.
- **Memory.** The input plus the k × M table of DP back-pointers, at most 8.4 MB for k ≤ 512.
- **Larger M.** At N = 1 million on continuous data, each doubling of M lowered the excess SSE
  about 3.9×, at little time cost. A larger M such as 8,192 is an option; the default stays 4096.

`lex_plain` (C++, Python, R; `Lex.lexPlain` in Scala) is the same method with the coarse start off.

## Deck preset (`lex_plain` with `deck()`)

This is the 2024 deck setting. Sorted values are grouped by a running-mean rule at precision 1e-4:
a value joins the current run if it is less than 1e-4 above the run's
mean (equal values always join). This **changes the
problem**: the groups are atoms, and the SSE is still reported on the raw values.

The preset also tries three starts (Collapse, equal-frequency, equal-width), keeps the one with the
lowest SSE, and refines it. It runs through `lex_plain`, so the coarse start is off: C++
`lex_plain(xs, k, deck())`, Scala `Lex.lexPlain(xs, k, Options.deck)`, Python and R
`lex_plain(x, k, precision = 1e-4, multi_start = TRUE)`. Its accuracy is limited by the grid.

## Tie rule (aligned with Ckmeans.1d.dp)

Among partitions with equal SSE, every solver returns the one whose cuts are **largest**, compared
from the last cut down. Ckmeans.1d.dp (R) does the same. How each part does it:

- **DP.** Scans split points right to left with a strict `<`, so the largest split index wins among
  exactly equal costs.
- **Refine, gain moves.** Among exactly equal best deltas, the move with the larger new bounds wins
  (right bound first).
- **Refine, tie sweeps.** After convergence, Refine applies moves that shift the window's bounds
  right with |δ| ≤ 4 ulp · mag, taking the largest new bounds, then runs gain sweeps again. It
  repeats this at most 8 rounds. Each tie sweep strictly raises the bound vector, so it terminates.
  With weights, only δ = 0 exactly counts as a tie.
- **`ckmeans_exact`.** Runs the same tie pass after the DP, so ties hidden by the DP's rounding also end
  on the largest cuts.
- **Collapse.** Still merges the leftmost pair on exactly equal cost. Refine's tie sweeps then
  move the result to the largest cuts.

Example: {0, 1, 2} with k = 2 gives sizes 2 | 1 (R: 2 | 1; Python `ckmeans`, which keeps the
smallest cut: 1 | 2).

## Weights (C++, Python and Scala)

Equal values are ordered by weight before the prefix sums, so results never depend on input order.
Atoms and groups carry double-double prefix sums of W = Σw, S = Σw·z and Q = Σw·z² (z the
centred value). Weighted counts replace counts in the merge cost, W₁W₂/(W₁+W₂)·(m₁ − m₂)², in the
shift delta and in the DP (bin cost Q − S²/W). A move whose piece takes the losing bin's whole
weight is costed from the bins' SSEs before and after. Zero weights are
allowed (the tests check that the exact SSE equals the problem with those points dropped);
negative, non-finite or all-zero weights are rejected. A bin whose values all have weight 0 has no
weighted mean: its reported centre is the weighted mean of all the data (centred value 0), and it
adds 0 to the SSE. Weighted `Atoms` do not build the unweighted prefix sums (`prefix_hi` /
`prefix_lo` are empty in C++, null in Scala).

Weighted SSE and centres are reported. Unit weights give the unweighted SSE up to rounding (tested).

**Extreme weight ratios.** The prefix sums are double-double, but the means, bin costs and move
deltas built from them are doubles. When one input mixes weights about 24 orders of magnitude
apart (the stress test uses 1e12, 1e-12 and 0), the small weights vanish from those quantities and
a move delta can have the wrong sign. Refine then flipped between two partitions forever. Each
Refine pass now also stops when a sweep ends on a partition already seen in that pass. This guard
never fires without a cycle: the 352 fixture runs are byte-identical with and without it. At such
ratios results are approximate. For weight ratios from 1e-6 to 1e6 `ckmeans_exact` matches a
brute-force optimum to 1e-7 relative (tested).

## Overflow (finite input too large for double arithmetic)

Finite values can still be too large for the sums LEX forms. Without a check such input gives a NaN
or infinite SSE and centres, or, worse, a finite but wrong partition. So every entry point
(`lex_hybrid`, `lex_plain`, `ckmeans_exact`, weighted or not, and R's `ckmeans()`) checks one rule
on the input before any solver runs. The solvers themselves (Atoms, Collapse, DP, Refine) have no
overflow checks, so input that passes gets exactly the results it got before the rule existed.

Let z = x − mean (the centred value the solvers use; with weights, the weighted mean),
Z = max |z|, DBL_MAX = 1.7976931348623157e308 and L = √DBL_MAX / 4 = 3.351951982485649e153.

- **Unweighted:** n · Z ≤ L.
- **Weighted:** max(W, 1) · Z ≤ L and W ≤ L, where W is the total weight.

A mean that is not finite (the sum of the values overflows, which needs values near ±1e308) makes
Z infinite or NaN, so it fails the rule as well. On failure all languages raise an error (C++
`std::invalid_argument`, Python `ValueError`, R `stop()`, Scala `IllegalArgumentException`) that
gives the limit L, the observed n · Z (or max(W, 1) · Z, or W), n or W, and Z, and says when the
mean is not finite. The check is one function per language (`check_range` in cpp/src/lex.cpp,
`checkRange` in Lex.scala; Python and R call the C++ one), run once after sorting and centring.

**Why the unweighted rule is enough (proof).** If n · Z ≤ L, then:
- the sum of squares SST = Σ z² ≤ n Z² ≤ (nZ)² ≤ DBL_MAX/16, and every bin SSE, DP cost and the
  reported SSE is at most SST;
- every run sum |S| = |Σ z| over a run of values is at most nZ ≤ L, so S² ≤ DBL_MAX/16 in the DP
  cost Q − S²/n;
- a merge cost n_a n_b/(n_a + n_b) (m_a − m_b)² is below n · 4Z² ≤ DBL_MAX/4;
- each Refine shift's mag (its rounding bound) is below 4 n² Z² (w₂ ≤ k, w₁ ≤ k n₁,
  k(n₁ + 1) < n², |g| + |m| ≤ 2Z), so a move's mag is below 8 (nZ)² ≤ DBL_MAX/2.
Z is the largest computed |z| (x − mean rounds monotonically in x), and every bound above has a
factor of at least 2 to spare, which covers the rounding of the sums.

**Weighted.** If max(W, 1) · Z ≤ L and W ≤ L, the same steps bound SST = Σ w z² ≤ W Z² ≤ L²
(for W < 1, W Z² < Z² ≤ L²), run sums |Σ w z| ≤ W Z ≤ L, the weight products W₁ W₂ and W_p W₁
(≤ W² ≤ L²), merge costs and DP costs. The W ≤ L part is needed because those weight products are
formed even when Z is tiny. Refine's mag is the one term the rule does not bound in general. Moving
a piece of weight W_p out of a bin of weight W₁ uses W_p W₁/(W₁ − W_p), which is large when the
weight left behind, r = W₁ − W_p, is small. The term is at most 4 W² Z²/r ≤ DBL_MAX/(4r), so mag
stays finite whenever r ≥ 1. That holds, and the proof goes through, when every positive weight is at
least 1 and the weight sums are exact (for example integer weights with W ≤ 2⁵³). With smaller
weights, mag can overflow near the limit. Such a move is never applied, because every comparison
with an infinite mag is false. **Conjecture, not proved:** for weighted input that passes the rule,
an overflowing mag never changes the result. The probes below found no counterexample.

**The rule is conservative.** It is sufficient, not necessary, so some input above the limit would
have fitted correctly and is rejected anyway. For example {−a, a} with a = √(DBL_MAX/2) has a finite
SST, but n · Z = √2 · √DBL_MAX > L. Real data are nowhere near the limit: it needs values around 1e150.

**Measured (box, x86-64, C++).** Every accepted run below gave the same bounds, moves and sweeps as
the same data scaled by 2⁻⁶⁰⁰. Scaling by a power of 2 is exact, and the scaled data are far from
overflow.
- 12,000 runs of clusters at ±2^490 to ±2^515, a third of them weighted (weights below 2¹⁰, down
  to 2⁻¹⁰ and smaller; hybrid, plain, exact, deck): 5,904 accepted. Of the 6,096 rejected, 2,960 have a true SST above DBL_MAX.
- 12,000 unweighted runs of heavy duplicates far from the centre (the shapes where Refine's mag is
  largest; hybrid, plain, exact): 774 accepted.
- 60,000 weighted runs near the limit, a quarter of the weights between 2⁻⁵⁹ and 1 and the rest 1
  to 3, so moves can leave little weight behind (hybrid, plain, exact): 55,608 accepted.
The 3000-case fuzz and all fixture runs are unchanged.

## Complexity summary

| Mode | Time after the O(n log n) sort | Memory |
|---|---|---|
| `lex_plain` | O(d log d) Collapse + sweeps × O(k) | O(n) |
| `lex_hybrid` (coarse M = 4096) | O(d log d) + O(k M log M) + refine | O(n + k M) |
| deck | as plain, on the grid groups | O(n) |
| `ckmeans_exact` | O(k d log d) | O(n + k d) |

## Implementations and intended differences

- **C++** (`cpp/`) is the production core. The Python bindings call it directly.
- **Scala** (`scala/`) is the reference implementation of the same algorithm. On the shared fixture
  set it is bit-identical to C++ for `lex_hybrid`, `lex_plain`, deck and `ckmeans_exact`, weighted and
  unweighted: cuts, SSE and centre bit patterns, start, atom, sweep and move counts.
  `CrossCheckTest` checks this on 356 runs: 64 weighted (some with zero weights, and one 17,000-value
  set that takes the coarse path over groups), 16 with k clamped, and 4 on `fma_60`, a case whose
  SSE changes if the C++ core is built with FMA contraction (x86-64, `-ffp-contract=fast
  -march=haswell`, at -O2 and -O3; the other cases miss the -O3 build). The Scala code mirrors the
  C++ files one to one (each Scala file names its C++ counterpart), including the 4-ary heap, the
  compact weighted prefixes and Refine's per-window reuse of the first shift.

Intended differences:
- Scala has no presorted-input flag (`Arrays.sort` always runs). It names the result fields
  `counts` / `means` (C++: `sizes` / `centers`), and its error messages spell doubles the Java way
  (`1.0E-4`).
- C++ is built without FMA contraction (`-ffp-contract=off` on GCC/Clang, `/fp:precise` on MSVC 2022+,
  `/fp:strict` on older MSVC; `cpp/cmake/LexFloatingPoint.cmake`). The JVM never contracts. This is
  what makes bit-identity possible on arm64.
