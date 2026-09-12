# IGNITE

**Immune Genomics and fuNctional Integration for Target Enrichment**

Code and data accompanying the manuscript *Genome-scale perturbation signatures from primary human CD4+ T cells improve genetics-based prioritization of immune drug targets*.

A positive-unlabeled (PU) learning model ranks immune drug-target candidates across **19,502 genes** and **24 features**, integrating rare- and common-variant genetic evidence and evolutionary constraint with perturbational functional genomics (genome-scale CRISPRi perturb-seq in primary human CD4+ T cells) and observational T-helper differential expression.

---

## Quickstart

```bash
make setup            # provision both stacks into .venv/ and .rlib/
make submission       # clone -> final_outputs/, offline, in one command
```

`make submission` runs the whole chain: version check, `tables`, `figures`, `verify`,
`final-outputs`. To work a stage at a time instead:

```bash
make check-versions   # check Python and R against the pins
make test             # fast contract checks, no fitting or rendering
make figures          # ~15 s: render every display item from committed tables
make verify           # did it reproduce?
make final-outputs    # assemble the submission package
```

`make figures` works on a **fresh clone with no model fitting**, since the table layer is committed. To rebuild that layer: `make tables` (~13 min), then `make verify-tables`.

`make help` lists everything. Interpreters are overridable: `make tables PY=... RSCRIPT=...`.

**Requirements.** Python 3.12 or later, with the pins in `requirements.txt`; R 4.4.1, with the pins in `R-requirements.txt`. `make setup` provisions both. The build runs on macOS and Linux (including WSL on Windows). macOS is needed only for the PDF versions of the figures, which R writes through macOS's quartz device; on Linux the build skips the PDFs and reproduces every table and number, while PNG pixels may differ from the macOS baseline because fonts and rasterization libraries differ. See `REPRODUCIBILITY.md`, PDF outputs and the supplementary target. On Debian/Ubuntu, `make setup` needs `build-essential` and the headers for freetype, harfbuzz, fribidi, libpng, libjpeg, libtiff, fontconfig, libcurl and libxml2; it names the exact `apt-get` command if any are missing.

---

## Two stages, with a committed contract

```
data/          committed inputs
  └─ stage 1: src/analysis/*.py  ──►  figure_data/   (COMMITTED — the contract)
                                        └─ stage 2: src/figures/*.R ──► final_plots/
```

| | |
|---|---|
| **Stage 1** | 30 Python steps: fit the models, run the tests, write 92 tables into `figure_data/`. **Step order is load-bearing — do not run with `-j`.** The `tables:` target comment lists the constraints. |
| **Stage 2** | 10 R renderers reading `figure_data/` and *nothing else*. **Renderers are pure** — every count, star, p-value and threshold is precomputed by stage 1. A figure needing a new number needs a stage-1 change. |

`figure_data/` is committed, which makes stage 2 standalone. As a result, `make tables` leaves the tree dirty under `figure_data/`, and `make verify-tables` is the check that it reproduced: a clean tree means the numbers came back byte-identical. `figure_data/PROVENANCE.md` maps every file to the script that writes it.

---

## What reproduces, and what does not

All results reported in the manuscript and supplementary tables were generated on macOS (Apple Silicon) using the environment provided in this repository. Minor numerical differences may occur when the analyses are reproduced on other operating systems or hardware, because of differences in numerical libraries, floating-point operations and parallel computation. These differences do not affect the qualitative conclusions or statistical interpretation of the analyses. In the exploratory discordance analysis (Supplementary Fig. S2 and Supplementary Table 9), the core gene set and features near the FDR threshold may differ. The measured differences on Linux are in `REPRODUCIBILITY.md`, Cross-platform reproduction.

| | Claim | Behaviour |
|---|---|---|
| `make verify-tables` | Did the **numbers** reproduce? | **Hard failure.** Portable across machines. The check to trust. |
| `make verify-figures` | Did the **pixels** reproduce? | Reported; never fails for rendered PNGs. Fatal only for absent files. |

`verify-figures` is soft because `ragg` rasterizes glyphs through freetype/libpng, which neither CRAN nor a lockfile can pin. **A different machine can reproduce every number exactly and still fail the PNG hashes.**

**Three things a green gate does not cover:**

- **The permutation null.** `make tables` does not rewrite the `label_permutation_*` files, so `verify-tables` confirms only that they are unchanged. They *were* regenerated under the current pins on 2026-09-06 against the 24-feature matrix — 1,000 permutations plus 99 scrambled-feature draws — so the committed draws are a product of this toolchain. What is still uncovered is repetition: they have been generated once, not shown to re-derive. `make permutation-null` (~8.5 h at `PERM_JOBS=11`, measured) regenerates them.
- **The opt-in targets.** `make tables` rewrites 81 of the 92 tables. The remaining 11 belong to targets that need network access or a large input (`permutation-null`, `discordance-network`, `vignette-panelc`, `gps-inputs`, `knn-signatures`) and pass `verify-tables` without being rewritten.
- **Rebuilding `final_plots/` does not test `figure_data/`.** 41 of the 92 tables are read by no renderer; most are sheets of the numbered supplementary tables. `verify-tables` is the check that matters; the figures follow from it.

The xgboost version pin fixes gene rank order: 3.2.0 versus 3.3.0 moves one gene from rank 9 to 11 deterministically, though aggregate metrics are unchanged. `_version_guard.py` runs at import as step 1 of `make tables` and catches a mismatch that would otherwise pass silently. `IGNITE_ALLOW_VERSION_MISMATCH=1` downgrades the check to a warning; outputs from such a run should not be committed.

---

## Reading the figures

### Sample sizes differ across panels — this is expected

A single drug-status rule (`furthest_stage`: Approved, In-trials for Phase 1/2/3, non-target for everything else) applies over each assay's set of measured genes. A gene with no cytokine reading does not appear in the cytokine panel. The caption clause "n = genes with data for this assay" documents this.

| Figure / panel | Approved | In-trials | Other | Denominator |
|---|---|---|---|---|
| Fig 2 **a** — IEI enrichment | 340 | 734 | 18,435 | 19,502-gene matrix (background only) |
| Fig 2 **b** — GWAS score | 340 | 727 | 18,435 | 19,502-gene matrix |
| Fig 2 **c** — missense z | 340 | 727 | 18,360 | matrix ∩ `mis.z_score` present |
| Fig 3 **b** — regulator residual | 206 | 318 | 9,744 | matrix ∩ regulator residual present |
| Fig 3 **c/d** — cytokine (+receptor) | 200 | 383 | 10,902 | cytokine-assay genes ∩ the universe |
| Th-subsets violin | 261 | 410 | 10,719 | genes with ≥1 Th-vs-Th0 DE value |
| Fig 4 **d** / score-by-group | 340 | 727 | 18,435 | 19,502-gene matrix |

- 734 versus 727: 734 counts every Phase 1/2/3 gene in the drug table minus the approved list, without universe intersection; 727 is the subset present in the feature matrix. The gap of seven genes is mitochondrially encoded metformin targets (MT-ND1 through MT-ND6), absent from the nuclear matrix. Panel a's background is the same 19,502-gene matrix; only its bar counts remain unintersected.
- Multi-gene drug targets are semicolon-joined (ruxolitinib maps to `JAK1;JAK2;JAK3;TYK2`). Producing scripts split on `;` and explode before the `groupby`; omitting that split undercounts the immune-indication counts.

## Display items

Seven main data figures (Fig 2-8) plus supplementary S1-S3. `tools/supplementary_figures.tsv` is authoritative for the S-numbers.

```bash
make figures       # ~15 s — render every display item from the committed tables
make final-outputs # assemble final_outputs/, numbered as the manuscript numbers them
```

| `final_outputs/` | Contents |
|---|---|
| `main_figures/` | Fig 2–8 as PNG, plus a `pdf/` vector copy of each |
| `supplementary_figures/` | S1–S3 as PNG, plus a `pdf/` vector copy of each |
| `supplementary_tables/` | the 15 numbered tables (`.csv` / `.xlsx`), straight out of `figure_data/` |
| `README.md` | generated manifest with a SHA-256 for every file |

`final_outputs/` is gitignored: it is assembled from committed content on demand, so committing it would create a second copy that can drift from the source.

### Rebuilding only the supplementary set

`make final-outputs` depends on `verify`, so it does not run until all ten display items are present and passing, which suits a submission bundle. Iterating on a single supplementary figure requires a full `make figures` and a rebuild of the main figures to see a changed S2.

```bash
make supplementary # 5 renderers; writes supplementary/; main figures untouched
```

It packages the same 15 tables and the same S1-S3 images into a separate tree because `final-outputs` opens with `rm -rf final_outputs`. A supplementary set assembled into that tree would be removed on the next full build without a record of its existence, and a partial set left in the submission folder would appear current.

**The two trees are not interchangeable.** They name the figures differently: `Supp Figure 1 - <title>.png` here versus `S1_<slug>.png` in `final_outputs/`, and neither carries both forms. Both are gitignored, and `make clean` removes both.

**What reproduces byte-for-byte.** On a matched toolchain, every `.png` and `.csv` in the package is byte-identical across rebuilds. The `.xlsx` tables are content-identical: every workbook part hashes the same except `docProps/core.xml`, which carries a wall-clock stamp. The `.pdf` copies differ only in `/CreationDate`, `/ModDate` and the `/ID` derived from them; `verify-figures` classifies PDFs as soft on that basis. `tools/compare_submission.py` compares two assembled packages under these rules.

---

## Data availability

Every input `make tables` needs is committed under `data/`; the offline build works from a clean clone. Per-file provenance is in `figure_data/PROVENANCE.md`.

**One exception, by size.** The perturb-seq knockdown signature matrices (`data/perturbseq/knn/signatures_*.npy`, 294 MB for the one that is read) are not committed and are regenerated instead:

```bash
make knn-signatures   # OPT-IN, network, ~3-4 min
```

which streams them from the public source via HTTP range requests. The single step that reads them, `src/analysis/build_vignette_panelD_knn.py`, backing panel d of the STAT4 vignette, skips with a message when they are absent and leaves its two committed `figure_data/` tables in place, so a clean clone still runs `make tables` end to end and passes `make verify-tables`. The row and column gene-name lists beside the matrices are committed, so their identity can be checked without them.

Upstream source: the genome-scale CRISPRi Perturb-seq resource in primary human CD4+ T cells of Zhu, Dann et al. (2025), released via the CZI Virtual Cells Platform (`s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`). That resource is the upstream scientific data source for every perturb-seq feature here and must be cited in any work using them.

A handful of files have no producing script and cannot be rebuilt from this repository: the external lists in `data/prospective/`, two files in `data/specificity/`, and the hand-curated 16-phecode list in `data/reference/immune_phecodes_curated.csv`. `figure_data/PROVENANCE.md` records which files these are.

### Third-party data

Redistributed here under the terms of their respective sources, each to be cited independently of this repository: **HGNC** (gene symbol reference), **Open Targets** (drug and association evidence; the live-query target is release 26.06), **GWAS Catalog**, **STRING** and **Gene Ontology** (S2's network and enrichment caches), and **gnomAD** (constraint metrics).


---

## Repository layout

| Path | |
|---|---|
| `data/` | committed inputs |
| `src/analysis/` | stage 1 — 30 Python steps. `README.md` there is the PU-model settings census |
| `src/figures/` | stage 2 — 10 renderers. `palette.R` is the shared colour convention |
| `figure_data/` | the committed contract. `PROVENANCE.md` maps every file to its producer |
| `outputs/` | stage-1 scratch no figure reads; nothing here is committed |
| `tools/` | environment setup, verification, packaging |
| `tests/` | pre-flight contract checks, run by `make test` |
| `final_plots/` | stage-2 output (mostly gitignored — rebuild with `make figures`) |
| `final_outputs/` | assembled on demand by `make final-outputs` |

`REPRODUCIBILITY.md` documents what the build reproduces, what it does not, and every
constraint affecting whether a re-run yields the committed numbers.
`DATA_AVAILABILITY.md` lists every external data source, what is committed, what is not,
and the terms each source is redistributed under.

## Citation

Cite the manuscript this repository accompanies, the upstream perturb-seq resource (Zhu, Dann et al., 2025, released via the CZI Virtual Cells Platform), and the third-party data sources named above, each under its own terms.

## License

MIT, see `LICENSE`. The license covers the code in this repository; redistributed third-party data remains under its own terms. Upstream data sources and attributions are in `NOTICE.md`.
