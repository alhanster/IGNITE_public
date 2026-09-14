# IGNITE

**Immune Genomics and fuNctional Integration for Target Enrichment**

Code and data accompanying the manuscript *Genome-scale perturbation signatures from primary human CD4+ T cells improve genetics-based prioritization of immune drug targets*.

A positive-unlabeled (PU) learning model ranks immune drug-target candidates across 19,502 genes and 24 features, integrating rare- and common-variant genetic evidence and evolutionary constraint with genome-scale CRISPRi perturb-seq in primary human CD4+ T cells and observational T-helper differential expression.

# Requirements

- Python >= 3.12, packages pinned in `requirements.txt`
- R 4.4.1, packages pinned in `src/figures/R-requirements.txt`
- macOS or Linux, including WSL on Windows. Vector PDF output requires macOS (quartz device); other platforms write PNG only.
- On Debian/Ubuntu: `build-essential` and headers for freetype, harfbuzz, fribidi, libpng, libjpeg, libtiff, fontconfig, libcurl, libxml2. `make setup` names the exact `apt-get` command if any are missing.
- Reference platform is Apple Silicon macOS. The build completes on Linux x86_64 and reproduces every table and number there, with model-derived values differing in the last digits; see `REPRODUCIBILITY.md`.

# Usage

```
make setup            # provision both stacks into .venv/ and .rlib/
make submission       # version check -> tables -> figures -> verify -> final-outputs
```

One stage at a time:

| Target | Time | Does |
|---|---|---|
| `make check-versions` | seconds | checks Python and R against the pins |
| `make test` | seconds | contract checks; no fitting or rendering |
| `make tables` | ~13 min | stage 1: 30 Python steps write `figure_data/` |
| `make figures` | ~15 s | stage 2: 10 R renderers write `final_plots/` |
| `make verify` | seconds | `verify-tables` + `verify-figures` |
| `make final-outputs` | seconds | assembles `final_outputs/` |
| `make supplementary` | seconds | renders and packages the supplementary items only |
| `make clean` | seconds | removes build outputs |
| `make help` | | lists every target |

- `make figures` runs on a fresh clone with no model fitting, because `figure_data/` is committed.
- Interpreters are overridable: `make tables PY=... RSCRIPT=...`
- Opt-in targets (network or long compute) are listed in `REPRODUCIBILITY.md`.

# Pipeline

```
data/          committed inputs
  └─ stage 1: src/analysis/*/*.py ─►  figure_data/   (committed)
                                        └─ stage 2: src/figures/*.R ──► final_plots/
```

- **Stage 1** — 30 Python steps fit the models, run the tests, and write 98 tables into `figure_data/`. Step order is load-bearing; do not run with `-j`.
- **Stage 2** — 10 R renderers read `figure_data/` and nothing else. Every count, star, p-value and threshold is precomputed by stage 1.
- `make tables` rewrites `figure_data/` in place, so the working tree is dirty after a run. `make verify-tables` is the check that it reproduced.
- `figure_data/PROVENANCE.md` maps every table to the script that writes it.

# Repository layout

| Path | Contents |
|---|---|
| `data/` | committed inputs |
| `src/analysis/` | stage 1 — 30 Python steps, one subfolder per display item plus `common/` helpers |
| `src/figures/` | stage 2 — 10 renderers plus `palette.R` and `vector_output.R` |
| `figure_data/` | the committed stage-1 output read by stage 2 |
| `outputs/` | stage-1 scratch; not committed |
| `tools/` | environment setup, verification, packaging, contract tests |
| `final_plots/` | stage-2 output; rebuild with `make figures` |
| `final_outputs/` | submission package; assembled by `make final-outputs` |

Display items: main Figures 2-8 and Supplementary Figures S1-S3, plus 15 numbered supplementary tables. `tools/supplementary_figures.tsv` and `tools/supplementary_tables.tsv` are authoritative for the numbering.

## tools/

| File | Does |
|---|---|
| `setup_env.sh` | provisions both stacks into `.venv/` and `.rlib/` (`make setup`) |
| `detect_python.sh` | selects an interpreter, preferring one with the pinned xgboost installed |
| `check_r_versions.R` | checks the R version and package pins (`make check-versions`) |
| `test_reproduction_contracts.py` | the contract checks run by `make test` |
| `run_step.sh` | runs one stage-1 or stage-2 step exactly as the build does |
| `write_baseline.sh` | regenerates `tools/final_plots.sha256` from a clean clone |
| `package_supplementary_tables.py` | packages the numbered tables from `figure_data/` |
| `package_supplementary_figures.py` | packages the supplementary figures |
| `write_submission_readme.py` | generates `final_outputs/README.md` with a checksum per file |
| `compare_submission.py` | compares two assembled `final_outputs/` packages |
| `supplementary_tables.tsv`, `supplementary_figures.tsv` | the numbering manifests |
| `final_plots.sha256` | the figure checksum baseline checked by `make verify-figures` |

# Data

Every input `make tables` needs is committed under `data/`; the offline build works from a clean clone. Per-file provenance for the inputs is in `data/DATA_AVAILABILITY.md`; for the derived tables, in `figure_data/PROVENANCE.md`.

One exception by size: the perturb-seq knockdown signature matrices (`data/perturbseq/knn/signatures_*.npy`) are not committed. Regenerate with `make knn-signatures` (network, ~3-4 min).

Third-party data redistributed here, each to be cited independently: **Open Targets** (release 26.06), **HGNC**, **GWAS Catalog**, **STRING** (v12.0), **Gene Ontology / UniProt-GOA**, **gnomAD** (v4.1.0), **IUIS** IEI gene list, the **GPS** dataset, **ChEMBL** and **ClinicalTrials.gov**, and the perturb-seq resource of Zhu, Dann et al. (2025). See `data/DATA_AVAILABILITY.md`.

# Documentation

| File | Contents |
|---|---|
| `REPRODUCIBILITY.md` | version pins, stage-1 step order, verification targets, opt-in targets, model configuration |
| `data/DATA_AVAILABILITY.md` | every external input, where it lives, what is committed, redistribution terms |
| `figure_data/PROVENANCE.md` | one line per committed table: the script that writes it |
| `src/analysis/README.md` | stage-1 layout and the scripts that train the PU model |
| `src/figures/README.md` | renderer-to-figure map and shared helpers |

# Citation

Cite the manuscript this repository accompanies, the upstream perturb-seq resource (Zhu, Dann et al., 2025, released via the CZI Virtual Cells Platform), and the third-party data sources named in `data/DATA_AVAILABILITY.md`, each under its own terms.

# License

MIT; see `LICENSE`. The license covers the code in this repository. Redistributed third-party data remains under its own terms.
