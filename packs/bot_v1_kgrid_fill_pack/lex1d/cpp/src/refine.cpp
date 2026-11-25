// 3-bin-window refinement. For the window (A, X, B) around bin X, each of X's two bounds can stay,
// give X's edge piece to the neighbour, or take the neighbour's edge piece: 3 x 3 - 1 = 8 moves.
// A piece is the part of one atom inside one bin, so moves never split an atom further.
//
// A sweep visits windows from right (X = k-2) to left (X = 1) and applies each window's best move
// immediately if it lowers SSE; among exactly equal deltas the move with the larger new bounds wins
// (compared right bound first). Sweeps repeat until none applies. With k = 2 the single window is
// (A = 0, X = 1). Bounds are updated in place.
//
// Ties (as Ckmeans.1d.dp: largest cuts). After convergence, tie sweeps apply moves that shift the
// window's bounds right (larger, right bound first) and whose |delta| <= kTieEps * mag (equal SSE up
// to rounding; exactly 0 with weights), taking the largest new bounds. If any applied, gain sweeps
// run again, and so on (at most kTieRounds rounds). Tie sweeps are not counted in RefineStats::sweeps; their moves are.
//
// Termination: a gain move is applied only if its computed delta is below -tolerance * mag, where mag
// bounds the size of the delta's terms. Rounding error is far below tolerance * mag, so every applied
// gain move strictly lowers the true SSE; no partition can repeat, and there are finitely many. Each
// tie sweep run strictly raises the bound vector (compared from the last bound), so it ends too.
// The rounding argument fails when weights span more than double's ~16 digits (say 1e12 and 1e-12 in
// one row): range sums then lose the small weights and gains can come out with the wrong sign. As a
// guard, a pass also stops when a sweep ends on a partition already seen in that pass. Without a
// cycle the guard never fires, so it changes no result.
#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <set>

namespace lex::detail {
namespace {

//                         move: 0  1  2  3  4  5  6  7
constexpr int kLeftOp[8]  = {1, 1, 0, 2, 0, 2, 1, 2}; // 0 stay, 1 give, 2 take
constexpr int kRightOp[8] = {1, 0, 1, 0, 2, 1, 2, 2};

struct Refiner {
    const Atoms& a;
    std::vector<int>& b;
    double tolerance;
    double tie_eps = a.weighted ? 0.0 : kTieEps; // weighted: only exactly zero deltas are ties
    double mag = 0; // accumulated by shift()

    /// SSE change from moving the bound between bins lo..mid-1 and mid..hi-1 to `to` (one piece of
    /// k elements with value g). Closed form: k n2/(n2+k) (g-m2)^2 - k n1/(n1-k) (g-m1)^2,
    /// where 1 is the bin losing the piece and 2 the bin gaining it.
    double shift(int lo, int mid, int to, int hi) {
        if (a.weighted) {
            if (to == mid) return 0.0;
            const int p0 = std::min(to, mid), p1 = std::max(to, mid);
            const double Wp = a.wsum(p0, p1);
            if (Wp == 0) return 0.0;
            const bool right = to > mid;
            const double W1 = right ? a.wsum(mid, hi) : a.wsum(lo, mid);
            const double W2 = right ? a.wsum(lo, mid) : a.wsum(mid, hi);
            const double m1 = right ? a.mean(mid, hi) : a.mean(lo, mid);
            const double m2 = right ? a.mean(lo, mid) : a.mean(mid, hi);
            if (W1 > Wp) {
                const double g = a.value[a.atom_id[p0]];
                const double ww1 = Wp * W1 / (W1 - Wp);
                const double ww2 = Wp * W2 / (W2 + Wp);
                const double s1 = std::abs(g) + std::abs(m1);
                const double s2 = std::abs(g) + std::abs(m2);
                mag += ww2 * s2 * s2 + ww1 * s1 * s1;
                return ww2 * (g - m2) * (g - m2) - ww1 * (g - m1) * (g - m1);
            }
            const double before = a.range_sse(lo, mid) + a.range_sse(mid, hi);
            const double after = a.range_sse(lo, to) + a.range_sse(to, hi);
            mag += std::abs(before) + std::abs(after);
            return after - before;
        }
        if (to == mid) return 0.0;
        const double k = std::abs(to - mid);
        const double g = a.value[a.atom_id[std::min(to, mid)]];
        const bool right = to > mid; // piece moves from the right bin to the left bin
        const int n1 = right ? hi - mid : mid - lo; // bin losing the piece
        const double m1 = right ? a.mean(mid, hi) : a.mean(lo, mid);
        const int n2 = right ? mid - lo : hi - mid; // bin gaining it
        const double m2 = right ? a.mean(lo, mid) : a.mean(mid, hi);
        const double w1 = k * n1 / (n1 - k);
        const double w2 = k * n2 / (n2 + k);
        const double s1 = std::abs(g) + std::abs(m1);
        const double s2 = std::abs(g) + std::abs(m2);
        mag += w2 * s2 * s2 + w1 * s1 * s1;
        return w2 * (g - m2) * (g - m2) - w1 * (g - m1) * (g - m1);
    }

    /// Gain mode: applies the best move around bin x if it lowers SSE by more than tolerance * mag.
    /// Tie mode: among moves that shift the bounds right with |delta| <= tie_eps * mag (tie_eps is
    /// kTieEps unweighted, 0 weighted; `tolerance` is not used here), applies the one with the
    /// largest new bounds.
    bool window(int x, bool tie) {
        const bool has_b = x + 2 < static_cast<int>(b.size());
        const int al = b[x - 1], l = b[x], r = b[x + 1];
        const int br = has_b ? b[x + 2] : r;
        const auto& id = a.atom_id;
        const int x_pieces = id[r - 1] - id[l] + 1;
        const int a_pieces = id[l - 1] - id[al] + 1;
        const int b_pieces = has_b ? id[br - 1] - id[r] + 1 : 0;
        // new position of each bound per op: stay, give, take
        const int give_l = std::min(a.start[id[l] + 1], r);
        const int take_l = std::max(a.start[id[l - 1]], al);
        const int give_r = std::max(a.start[id[r - 1]], l);
        const int take_r = has_b ? std::min(a.start[id[r] + 1], br) : r;

        // The first shift of the usual order depends only on the left op, so its delta and its mag
        // contribution are computed once per left op (mag starts at 0 for each move, so 0 + m1 == m1
        // exactly and results are unchanged). The give-left/take-right move shifts in the other order.
        double first_delta[3], first_mag[3];
        bool first_done[3] = {false, false, false};
        int best = -1;
        double best_delta = 0, best_mag = 0;
        int best_l = 0, best_r = 0;
        for (int j = 0; j < 8; ++j) {
            const int lo = kLeftOp[j], ro = kRightOp[j];
            const int gives = (lo == 1) + (ro == 1);
            const bool viable = (lo != 2 || a_pieces >= 2) && (ro == 0 || has_b) && (ro != 2 || b_pieces >= 2) &&
                                (lo == 2 || ro == 2 || x_pieces > gives);
            if (!viable) continue;
            const int nl = lo == 0 ? l : lo == 1 ? give_l : take_l;
            const int nr = ro == 0 ? r : ro == 1 ? give_r : take_r;
            mag = 0;
            // apply X's "take" first so X is never empty mid-move
            double delta;
            if (lo == 1 && ro == 2) {
                delta = shift(l, r, nr, br);
                delta += shift(al, l, nl, nr);
            } else {
                if (!first_done[lo]) {
                    first_delta[lo] = shift(al, l, nl, r);
                    first_mag[lo] = mag;
                    first_done[lo] = true;
                }
                mag = first_mag[lo];
                delta = first_delta[lo];
                delta += shift(nl, r, nr, br);
            }
            const bool larger = best >= 0 && (nr > best_r || (nr == best_r && nl > best_l));
            if (tie) {
                if (!(nr > r || (nr == r && nl > l)) || !(std::abs(delta) <= tie_eps * mag)) continue;
                if (best < 0 || larger) { best = j; best_l = nl; best_r = nr; }
                continue;
            }
            // equal deltas: the move whose new bounds are larger, right bound first (exact equality)
            if (best < 0 || delta < best_delta || (delta == best_delta && larger)) {
                best = j; best_delta = delta; best_mag = mag; best_l = nl; best_r = nr;
            }
        }
        if (best >= 0 && (tie || best_delta < -tolerance * best_mag)) {
            b[x] = best_l;
            b[x + 1] = best_r;
            return true;
        }
        return false;
    }
};

} // namespace

RefineStats refine(const Atoms& a, std::vector<int>& b, double tolerance,
                   const std::function<void(const std::vector<int>&)>& on_move) {
    RefineStats stats;
    const int k = static_cast<int>(b.size()) - 1;
    if (k < 2) return stats;
    Refiner r{a, b, tolerance};
    const int top = k == 2 ? 1 : k - 2;
    auto run = [&](bool tie) { // sweeps until none applies; true if any move was applied
        bool any = false;
        std::set<std::vector<int>> seen; // partitions at the end of each sweep (cycle guard)
        for (bool improved = true; improved;) {
            if (!tie) ++stats.sweeps;
            improved = false;
            for (int x = top; x >= 1; --x) {
                if (r.window(x, tie)) {
                    improved = any = true;
                    ++stats.moves;
                    if (on_move) on_move(b);
                }
            }
            if (improved && !seen.insert(b).second) break; // cycle: only with extreme weight ratios
        }
        return any;
    };
    run(false);
    for (int round = 0; round < kTieRounds && run(true); ++round)
        if (!run(false)) break;
    return stats;
}

} // namespace lex::detail
