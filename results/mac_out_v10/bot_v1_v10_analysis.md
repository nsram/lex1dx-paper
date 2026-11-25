# v10 Mac pack analysis

Sriram ran `bot_v10_rerun_mac_pack` on Mina (macOS 15.8.1, arm64 T6031, Apple clang 17.0.0, `caffeinate -dims`, `OMP_NUM_THREADS=1`). Same machine as v8/v9. lex1d SHA `8a691ee005b657f985ffa299c57542b79570b08c` matches the freeze.

Scaling started 11:10:31 PM ET and finished around 11:28:37 PM ET (~18 min). Controls started 11:28:37 PM ET and finished 11:44:50 PM ET (16.2 min, same wall as v9). `all.sh` was about 34 min, well under the ~1 h + 16 min estimate. v8 scaling was also ~18 min; the hour figure was a conservative estimate, not the measured wall.

The three result zips agree where they overlap: scaling and controls `runs.csv` / `summary.md` are bit-identical across copies. The combined zip's `out/env.txt` is the controls overwrite (11:28:37 PM ET). The scaling-start env is in the scaling-only zip. Attachments were copied, not modified. Preferred tree: `mac_out_v10/out/` plus `_from_scaling_zip/out/env.txt`.

## Completeness

| part | cells | status ok | expected |
|---|---:|---:|---:|
| scaling (v10) | 270 | 270 | 270 |
| scaling (v8) | 270 | 270 | 270 |
| controls (v10) | 64 | 64 | 64 |
| controls (v9) | 64 | 64 | 64 |

Scaling is 27 datasets × (1 sort + 3 methods × 3 k). Controls are the v9 A/B/C grid (3 × 8M + audio; quantile at k=512). No crashes, no incomplete rows.

## k-fallback

It did not fire. All 81 scaling sig5 cells and all 18 controls sig5 cells have `sig=5`. Smallest D′ on the scaling grid is 48,537 (lognormal 2^16), far above k=512. These families never need a 6th digit or a raw-atom fallback.

## sig5 vs v8 (81 matched cells)

Raw Refine changes quality a little and cost a lot. Group, atoms, Collapse, and DP times are unchanged (median v10/v8 ≈ 1.00). Almost all of the extra time is Refine (median 8.5×, range 2.3–128×), because moves now walk raw distinct-value atoms of size n.

**REE.** 33 cells better, 48 worse, 0 exactly equal. ΔREE (v10−v8) median +2.9×10⁻⁵, range −0.00228 to +0.00116. No cell has v10 REE > 0.01%. Grid max is 0.00526% (normal 2^18, k=256). Largest improvement is uniform 8M k=512 (0.00358% → 0.00130%). Largest regression is uniform 2^19 k=256 (0.00216% → 0.00333%). This matches the earlier local A/B: raw Refine is not strictly never-worse, and the changes are tiny.

**Time.** post_s ratio v10/v8 median 1.23, range 0.96–2.63. The high end is the 8M cells (uniform k=64: 0.056 s → 0.147 s in scaling). Small-n cells are essentially unchanged.

**Peak RSS.** Δ median +44 MB, range −12 to +398 MB. v8 range 8.9–204 MB; v10 range 9.5–497 MB. The one decrease (uniform 2^17 k=512: 27.5 → 15.3 MB) is small-n peak-RSS noise. At 8M, sig5 jumps from ~99–137 MB to ~421–497 MB in scaling (controls B: 419–552 MB). That is the raw Atoms prefix sums of length n, on top of the grouped atoms.

**Refine moves.** Median 836 → 5,301 (ratio median 5.23, range 0.92–58). At 8M, moves go from ~1k–6k to ~32k–102k, and sweeps rise with them.

## Tuned and Grønlund (scaling)

Tuned SSE is bit-identical to v8 on 81/81 cells. Wilber (patched) is bit-identical on 81/81. Time ratios are median 1.000 (tuned, range 0.95–1.30) and 0.999 (wilber, range 0.92–1.34). Medians are noise. A few large cells were 10–34% slower (uniform 8M tuned 1.23×; normal 2^19 k=512 wilber 1.34×), in line with ordinary run-to-run scatter on this machine, not a code change.

## Controls vs v9

**A. Pristine vs patched.** Same picture as v9. Patched/pristine post_s is 0.991–1.035 (median 1.004). One SSE last-digit mismatch, same cell as v9: uniform 8M k=64, 170.65890856765512 vs 170.65890856764022 (relative 8.7×10⁻¹⁴). Every other A cell matches. v10 A vs v9 A: 24/24 SSE bit-identical, time ratio 0.970–1.013 (median 0.997).

**B. sig5exact.** Unchanged. 12/12 SSE bit-identical to v9, time ratio 0.988–1.007 (median 1.000). REE values match. Peak RSS jumps around (audio k=512 was 356 MB in v9 and 155 MB in v10; a few normal cells dropped ~70 MB). That is peak-RSS noise, not a pipeline change. sig5 itself is slower and fatter, as in the scaling grid (e.g. uniform 8M k=64: 0.056 s / 99 MB in v8 → 0.117 s / 419 MB in controls B).

**C. Quantile.** Tuned is bit-identical to v9 at all six (dataset, q) cells, times within 1%. sig5 REE is slightly better at all six (normal q=0.75: 0.00204% → 0.00132%; audio q=0.75: 0.00209% → 0.00191%) and 1.27–1.84× slower. Default q=0.75 is still a reasonable pick; the sensitivity pattern did not change.

Controls also reconfirm 14/14 tuned/wilber SSE bit-identical to v8.

## Paper headlines at 8M

Recomputed on the v10 controls session: sig5 and sig5exact from experiment B, pristine from experiment A, end-to-end = post_s + sort_s. Tuned times are from v10 scaling (tuned is not in A/B) against the same pristine clocks. The twelve cells are uniform / normal / lognormal 8M and audio, at k=64, 256, 512.

| method | post-sort vs pristine | end-to-end vs pristine | REE max | peak RSS |
|---|---|---|---|---|
| sig5 (v10) | **33.2–75.9×** | **8.5–17.3×** | **0.00371%** | **419–552 MB** |
| sig5 (v8/v9, same recipe) | 52.8–153.9× | 9.3–20.1× | 0.00453% | 99–137 MB |
| sig5exact (v10) | 9.4–32.5× | 5.3–13.4× | 0.00089% | 109–165 MB |
| sig5exact (v9) | 9.3–32.4× | 5.3–13.9× | 0.00089% | 109–356 MB (RSS noisy) |
| tuned (v10) | 16.4–32.9× | 6.9–14.2× | 0.00216% | 683–747 MB |
| tuned (v8/v9) | 15.9–33.3× | 6.8–14.1× | 0.00216% | 682–746 MB |

The old paper line (sig5 52–152× post-sort, 9.3–20× end-to-end) should be replaced with **33–76× post-sort and 8.5–17× end-to-end**. sig5exact (9.3–32×) and tuned (16–33×) can stay. REE is still far under 0.01%; the 8M max is a bit better than v8. RSS for sig5 is now in the mid-hundreds of MB at 8M, because Refine keeps raw prefix sums of length n.

The low end of the sig5 range is still normal 8M k=256, where pristine happened to be unusually fast (5.76 s). Same cell was the low end in v9 (52.8×).

Using scaling sig5 times against the same pristine clocks gives a similar range (33.0–69.4× post-sort, 8.5–17.3× end-to-end, RSS 421–497 MB). Prefer the controls B numbers: they sit in the same session as pristine, which is how the v9 headlines were built.

## Anomalies and caveats

- k-fallback never appeared on this grid. The two-cluster regression (and any case with D′ < k) is covered by tests, not by these Mac cells.
- uniform 8M k=64 sig5 post_s is 0.147 s in scaling and 0.117 s in controls B, a 25% same-night gap. Headlines use controls B.
- The controls `summary.md` is still titled “v9 controls” (pack leftover). The numbers in it are the v10 run.
- Combined `env.txt` overwrites the scaling-start env. Both are kept.
- One pristine-vs-patched SSE last-digit mismatch (uniform 8M k=64) is unchanged from v9 and is negligible.
- Peak RSS is jumpy on methods whose memory did not change (sig5exact, a couple of small sig5 cells). Treat RSS ranges as approximate.
- Tuned and wilber quality are unchanged. A handful of their times were 10–34% slower than v8; medians are 1.00.
- Raw Refine can be slightly worse or slightly better than group-edge Refine. Nothing on this grid exceeds 0.01% REE.
