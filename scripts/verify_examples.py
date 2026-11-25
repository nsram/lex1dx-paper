"""Exact-fraction checks for the ten-value pipeline example and the seven-number Refine example.
Does not read or trust older prose; all SSE values are recomputed."""
from fractions import Fraction as F
from itertools import combinations

def sse(parts):
    t = F(0)
    for b in parts:
        m = F(sum(b), len(b))
        t += sum((F(v) - m) ** 2 for v in b)
    return t

def ward(A, B):
    wa, wb = len(A), len(B)
    ma, mb = F(sum(A), wa), F(sum(B), wb)
    return F(wa * wb, wa + wb) * (ma - mb) ** 2

def collapse(atoms, target, tie="left"):
    G = [list(g) for g in atoms]
    log = []
    while len(G) > target:
        costs = [(ward(G[i], G[i + 1]), i) for i in range(len(G) - 1)]
        minc = min(c[0] for c in costs)
        cands = [c for c in costs if c[0] == minc]
        i = cands[0][1] if tie == "left" else cands[-1][1]
        log.append((minc, i, len(cands), G[i][:], G[i + 1][:]))
        G[i : i + 2] = [G[i] + G[i + 1]]
    return G, log

def eight_moves(parts):
    L, X, R = parts
    out = []
    for lc in (0, 1, 2):
        for rc in (0, 1, 2):
            if lc == 0 and rc == 0:
                continue
            L2, X2, R2 = list(L), list(X), list(R)
            ok = True
            if lc == 1:
                ok = bool(X2); L2, X2 = (L2 + [X2[0]], X2[1:]) if ok else (L2, X2)
            elif lc == 2:
                ok = bool(L2); X2, L2 = ([L2[-1]] + X2, L2[:-1]) if ok else (X2, L2)
            if ok and rc == 1:
                ok = bool(X2); R2, X2 = ([X2[-1]] + R2, X2[:-1]) if ok else (R2, X2)
            elif ok and rc == 2:
                ok = bool(R2); X2, R2 = (X2 + [R2[0]], R2[1:]) if ok else (X2, R2)
            if not ok or not (L2 and X2 and R2):
                out.append(((lc, rc), None, None))
            else:
                p = [L2, X2, R2]
                out.append(((lc, rc), p, sse(p)))
    return out

x = [9, 12, 14, 17, 19, 22, 26, 27, 32, 36]
atoms = [[v] for v in x]

opts = sorted((sse([x[:a], x[a:b], x[b:]]), a, b) for a, b in combinations(range(1, 10), 2))
assert opts[0][0] == 83 and opts[0][1:] == (4, 8)
assert sum(1 for t in opts if t[0] == 83) == 1

G5, log5 = collapse(atoms, 5, "left")
assert [list(g) for g in G5] == [[9, 12, 14], [17, 19], [22], [26, 27], [32, 36]]
# LEX heap prefers the leftmost pair on a cost tie (tuned_collapse.hpp: ka < kb || (ka == kb && a < b)).
assert log5[0][0] == F(1, 2)  # 26 with 27
assert all(step[0] == log5[i][0] for i, step in enumerate(log5))  # sanity

restricted = []
for a, b in combinations(range(1, 5), 2):
    parts = [sum(G5[:a], []), sum(G5[a:b], []), sum(G5[b:], [])]
    restricted.append((sse(parts), parts))
restricted.sort(key=lambda t: t[0])
assert restricted[0][0] == F(424, 5)  # 84.8
assert restricted[0][1] == [[9, 12, 14, 17, 19], [22, 26, 27], [32, 36]]
assert sum(1 for t in restricted if t[0] == F(424, 5)) == 1

dp = restricted[0][1]
moves = eight_moves(dp)
improving = [m for m in moves if m[2] is not None and m[2] < sse(dp)]
assert len(improving) == 1
assert improving[0][1] == [[9, 12, 14, 17], [19, 22, 26, 27], [32, 36]]
assert improving[0][2] == 83
assert sse([[9, 12, 14, 17, 19]]) - sse([[9, 12, 14, 17]]) == F(144, 5)  # 28.8
assert sse([[19, 22, 26, 27]]) - sse([[22, 26, 27]]) == 27
assert sse(dp) - 83 == F(9, 5)  # 1.8

G3, _ = collapse(atoms, 3, "left")
assert [list(g) for g in G3] == [[9, 12, 14], [17, 19, 22], [26, 27, 32, 36]]
assert sse(G3) == F(1081, 12)  # 90.0833...
G3r, _ = collapse(atoms, 3, "right")
assert [list(g) for g in G3r] == [list(g) for g in G3]  # same 3-bin end point
c3_moves = eight_moves(G3)
assert all(m[2] is not None and m[2] > sse(G3) for m in c3_moves)
assert abs(float(sse(G3) / 83 - 1) - 0.085) < 5e-4  # 8.5%

# seven-number paired-move example
p0 = [[3, 4], [6, 7, 9], [14, 22]]
p1 = [[3, 4, 6], [7, 9, 14], [22]]
assert sse(p0) == F(223, 6)
assert sse(p1) == F(92, 3)
assert sse(p0) - sse(p1) == F(13, 2)
sm = eight_moves(p0)
singles = [m for m in sm if (m[0][0] == 0 or m[0][1] == 0) and m[2] is not None]
assert all(m[2] > sse(p0) for m in singles)
paired_hit = [m for m in sm if m[1] == p1]
assert len(paired_hit) == 1 and paired_hit[0][2] == F(92, 3)

print("OK ten-value and seven-number arithmetic")
print("  optimum SSE 83 (unique)")
print("  DP on five groups SSE", float(F(424, 5)), "=", F(424, 5))
print("  Refine recovers 83; net change", float(F(-9, 5)))
print("  Collapse-to-3 SSE", float(F(1081, 12)), "=", F(1081, 12))
print("  all eight Refine moves from Collapse-to-3 raise SSE")
print("  seven-number 223/6 -> 92/3, drop 13/2")
