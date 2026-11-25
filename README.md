# lex1dx-paper

Replication material for the paper *LEX: Fast scalar quantization by grouping
and local exchanges* by N. Sriram (draft dated 25 November 2026, file prefix
`lex1dx_paper_2026-11-25_r5`, with supplement
`lex1dx_supplement_2026-11-25_r5`).

The package is archived at Zenodo, https://doi.org/10.5281/zenodo.24543220.

LEX is the C++ library [lex1dx](https://github.com/nsram/lex1dx), with a Scala
port and Python and R interfaces. This repository holds the Mac benchmark packs
that produced the measurements, the output CSVs of those runs, the script that
rebuilds every result figure and table from those CSVs, and the figures and
tables themselves. The paper TeX sources and PDFs are not included.

The paper itself has no run or version tags. Pack names, output folders, and
the cell-by-cell provenance below are the record of which run produced which
number.

## Layout

| Path | What it holds |
|---|---|
| `paper/figs/` | The eight result figures used in the main paper and the supplement. |
| `paper/tables/` | The eight generated table fragments used in the main paper and the supplement. |
| `scripts/v4_figs_tables.py` | Rebuilds those figures and tables from the packaged CSVs. |
| `scripts/verify_examples.py` | Rechecks the exact-fraction worked examples in the paper. |
| `run/` | Wrappers: `figs_tables.sh` and `verify_examples.sh`. |
| `results/` | Output of the Mac runs used in the paper, plus two unused or superseded runs kept for provenance (see below). |
| `packs/` | The Mac packs that produced those outputs. Each pack is self-contained: frozen lex1d, vendors, drivers, and `run/` scripts. Large generated inputs are not included; each pack rebuilds them. |
| `requirements.txt` | Python packages for `bash run/figs_tables.sh`. |

Tables 1 and 2 of the main paper (configurations and experimental design) are
written in the paper source, not by the script.

## Machine and environment

Every timing run in `results/` was done on one machine, Mina, an Apple Mac
(arm64, Darwin 24.6.0, macOS 15.8.1, build 24H32). Pack READMEs describe it as
an M3 Max. The environment recorded in each `out/env.txt` is:

- Python 3.13.18, NumPy 2.3.5
- Apple clang 17.0.0 (clang-1700.4.4.1), CMake 4.5.0
- R 4.6.2, Ckmeans.1d.dp 4.3.6 with `method="linear"` (Ckmeans packs only)
- `OMP_NUM_THREADS=1`, and `caffeinate` to keep the machine awake
- lex1d frozen at commit `8a691ee005b657f985ffa299c57542b79570b08c`
- exact solver: Grønlund et al.'s `kmeans1d` at `585386d`, in the patched
  `wilber-interp` form unless a row is labelled `wilber_pristine` (the
  unmodified vendor build)

Synthetic inputs were generated on that Mac by `gen_synth.py` with
`SeedSequence(2048, spawn_key=(stream_id, N))`. Hashes are in each pack's
`data/synth/MANIFEST.sha256`. The audio file is the same recording as in
lex1d-paper (`dharini-abhogi.txt.gz`, N = 8,168,256). GPT-2 tensors are fetched
by the scaling pack when needed. Protocol: one warm-up plus three timed
repetitions, median (upper-middle). Times are post-sort unless noted;
`std::sort` is logged separately as `sort_s`.

The default coarse-group budget is M = max(4096, 32k). It is 4096 for k ≤ 128,
so the k curves evaluate that default policy, not k alone.

## How to regenerate figures and tables

From the repository root, with Python 3.12+ and the packages in
`requirements.txt`:

```
bash run/figs_tables.sh
```

That overwrites `paper/figs/v4_*.pdf` and `paper/tables/v4_*.tex` from the
packaged CSVs. It does not read anything outside this repository. The
regenerated tables match the packaged fragments byte for byte. Seven of the
eight figures match after PDF date and ID fields are ignored. `v4_seven.pdf`
comes out slightly smaller under the script's current figsize (3.25 by
2.85 in) than the packaged file, which is the snapshot used in the paper.

To recheck the exact-fraction pipeline and Refine examples:

```
bash run/verify_examples.sh
```

## How to rerun a Mac pack

Do this on the same class of machine if you want comparable times. Unzip or
copy one directory from `packs/` to its own folder, add the audio file as
`data/dharini-abhogi.txt.gz` when the pack needs it, and from that pack folder
run the command in its README. The commands are:

```
bash run/scaling.sh          # packs/bot_v8_scaling_mac_pack
bash run/controls.sh         # packs/bot_v9_controls_mac_pack
bash run/all.sh              # packs/bot_v10_rerun_mac_pack (see note below)
bash run/ckmeans_matched.sh  # packs/bot_v13_ckmeans_matched_pack
bash run/kgrid_fill.sh       # packs/bot_v1_kgrid_fill_pack
bash run/kgrid_lex_fix.sh    # packs/bot_v2_kgrid_fill_pack
```

Do not pass paths. Each script finds `python3` and a C++ compiler itself. The
Ckmeans packs also find `Rscript` and install `Ckmeans.1d.dp` into a pack-local
`rlib/` if it is missing. Runs are resumable: cells already marked `ok` in the
output CSV are skipped. Do not run two timing packs at once on the same
machine, because that skews the times.

Outputs land in that pack's `out/` directory. The CSVs already in `results/`
are the ones the paper uses.

Supplement Table S7 variant timings were produced by `packs/gronlund_bench`
(`run_bench.py` and `driver.cpp`). Only the output used in the table is needed
to rebuild the figures and tables, and that output is in
`results/gronlund/mac_out_2026-11-24/`.

## Figure and table map

All result figures and tables are written by `scripts/v4_figs_tables.py`.
Schematic figures (`v4_pipeline.pdf`, `v4_seven.pdf`) use only the exact
fractions checked by `scripts/verify_examples.py`.

Main paper:

| Paper object | File | Source CSVs and rows |
|---|---|---|
| Figure 1 (pipeline) | `paper/figs/v4_pipeline.pdf` | Worked example, no CSV. |
| Figure 2 (Refine moves) | `paper/figs/v4_seven.pdf` | Worked example, no CSV. |
| Figure 3 (max REE vs N) | `paper/figs/v4_ree_max.pdf` | v8 `tuned` and `sig5` at k ∈ {64,256,512}; precision-grouped fill-in at k ∈ {16,32,128} from kgrid_lex_fix, scored against kgrid_fill wilber SSE*. Synthetics only. |
| Figure 4 (time and RSS vs N at k=512) | `paper/figs/v4_time_rss_k512.pdf` | v8 `wilber`, `tuned`, `sig5`, and `sort`, all at k=512. |
| Table 3 (headline, k=512) | `paper/tables/v4_headline_k512.tex` | Upper block, four large inputs: exact reference and grouped exact from v9 `wilber_pristine` and `sig5exact`; full-precision LEX from v8 `tuned`; precision-grouped LEX from v9 experiment B `sig5`. Lower block, two GPT-2 matrices: v8 `wilber`, `tuned`, and `sig5`. Their baseline is the weighted table-free solver, marked † in the table. |
| Figure 5 (D₅/D vs N) | `paper/figs/v4_dprime.pdf` | v8 diagnostic columns, plus fill-in k from kgrid_lex_fix. |
| Figure 6 (time vs k at N≈8M) | `paper/figs/v4_time_vs_k.pdf` | Ckmeans: latest ok sorted row per (dataset, k), see Ckmeans note. LEX and exact: v8 at k ∈ {64,256,512}, kgrid_lex_fix / kgrid_fill wilber at k ∈ {16,32,128}. |

Online supplement:

| Paper object | File | Source CSVs and rows |
|---|---|---|
| Table S1 (REE by family) | `paper/tables/v4_ree_summary.tex` | v8 `tuned` (81 fits) and v8 `sig5` plus kgrid_lex_fix fill-in (102 fits). |
| Figure S1 (REE vs N by k) | `paper/figs/v4_ree_detail.pdf` | Same sources as Figure 3, per k. |
| Table S2 (headline, k=64,256,512) | `paper/tables/v4_headline_full.tex` | Same sources as the upper block of Table 3, three k. |
| Table S3 (Collapse quantile) | `paper/tables/v4_quantile.tex` | v9 experiment C, `tuned` and `sig5`, q ∈ {0.5, 0.75, 0.9}, Normal 2²³ and audio, k=512. |
| Table S4 (GPT-2) | `paper/tables/v4_gpt2.tex` | v8 `wilber`, `tuned`, and `sig5` on `h17_c_attn` and `h35_c_proj`. Baseline is the patched table-free solver. |
| Figure S2 (speedup vs N) | `paper/figs/v4_speedup_vs_N.pdf` | v8 plus kgrid fill-in, same split as Figure 3. |
| Table S5 (Ckmeans vs LEX vs exact, N≈8M) | `paper/tables/v4_ckmeans_8M.tex` | Same rows as Figure 6, four large inputs, k ∈ {16,32,64,128,256,512}. |
| Table S6 (Ckmeans at N=2²⁰) | `paper/tables/v4_ckmeans_1M.tex` | Same construction as Table S5, three synthetics at N=2²⁰, k ∈ {16,32,64,128,256,512}. |
| Table S7 (kmeans1d variants, audio) | `paper/tables/v4_variants_audio.tex` | `results/gronlund/mac_out_2026-11-24/out/results.csv`, in-process `driver_s` for `wilber-interp`, `hirc`, `wilber-binary`, and `monotone:compute` at k=64, 256, 512. These are not the main-grid post-sort medians. |

The small-k check at k=4 and k=8 on Normal, N=2²⁰, is stated in the text of the
results section (17 and 19 times faster after the sort, 6.9 and 8.0 times
including it). Those ratios come from `mac_out_kgrid_lex_fix` (sig5) and
`mac_out_kgrid_fill` (wilber and sort) and are not in a table fragment.

CSV locations under `results/`:

- v8: `mac_out_v8/out/scaling/runs.csv` (pack `bot_v8_scaling_mac_pack`, 24 November 2026, 20:42 ET)
- v9: `mac_out_v9/out/controls/runs.csv` (pack `bot_v9_controls_mac_pack`, 24 November 2026, 21:24 ET)
- kgrid_fill: `mac_out_kgrid_fill/out/lex/runs.csv` and `out/ckmeans/runs.csv` (pack `bot_v1_kgrid_fill_pack`, 25 November 2026, 11:28 ET)
- kgrid_lex_fix: `mac_out_kgrid_lex_fix/out/lex/runs.csv` (pack `bot_v2_kgrid_fill_pack`, 25 November 2026, 12:20 ET)
- Ckmeans folders, see next section
- gronlund: `gronlund/mac_out_2026-11-24/out/results.csv` (24 November 2026)

## Ckmeans provenance

`scripts/v4_figs_tables.py` concatenates four CSVs, keeps rows with
`method=Ckmeans.1d.dp`, `status=ok`, and `input=sorted`, and takes the last
remaining row per (dataset, k). In that concatenation order, the paper cells
are:

- k ∈ {16, 64, 128, 256, 512}: `results/ckmeans_matched_k64_256/` (copy of
  `mac_out_v13_k64_256`, pack `bot_v13_ckmeans_matched_pack`, 25 November 2026,
  10:15 ET). The pack README still mentions an earlier v11 pack name. The
  recorded env has the same driver, with a 100 GiB memory cap.
- k ∈ {32}, and the small-k check at k=4 and 8 on Normal, N=2²⁰:
  `mac_out_kgrid_fill/out/ckmeans/runs.csv`.

Two earlier Ckmeans folders are included for provenance but are not the last
row for any paper cell:

- `results/ckmeans_matched/` is a copy of `mac_out_v11_partial`. It has ok
  rows at k=16 and 128, and skipped 8M cells at k=512 (memory cap).
- `results/ckmeans_matched_k512/` is a copy of `mac_out_v12_k512`. It is the
  first run that finished the 8M, k=512 cells. Those cells were run again in
  v13 and the v13 rows are the ones in the tables.

## Runs included but not used in the figures and tables

- `results/mac_out_v10/` and `packs/bot_v10_rerun_mac_pack`: a full rerun of
  the v8 grid and v9 controls with raw-data Refine and a k-preserving digit
  fallback. The paper reports group-level Refine, so the published numbers
  come from v8 and v9, not v10. See `results/mac_out_v10/NOTE.txt`.
- LEX (`sig5`) rows in `results/mac_out_kgrid_fill/out/lex/runs.csv`: that
  first k-grid LEX run used the wrong Refine path. Use
  `results/mac_out_kgrid_lex_fix` instead. The wilber and Ckmeans cells in
  the same kgrid_fill output are valid and are used. See
  `results/mac_out_kgrid_fill/NOTE.txt`.

## License

The code (`scripts/`, `run/`, `packs/`) is under the MIT License (see
`LICENSE`).

The figures, tables, and results data are under Creative Commons
Attribution 4.0 International (see `LICENSE-CC-BY-4.0`), © 2026 Sriram
Natarajan.

The audio recording used by the packs, when you add it, is CC BY 4.0 by its
original authors (Carnatic Varnam Dataset, doi:10.34810/DATA2311), as in
lex1d-paper.
