# Reproducibility

Version pins, step order, verification targets, and model configuration. Companion to the Makefile: the Makefile says how to run each stage, this says what each stage needs.

# Environment

`make setup` provisions both stacks into project-local, gitignored directories, `.venv/` and `.rlib/`. `make check-versions` validates both against the pins.

## Python — stage 1, `make tables`

Requires Python >= 3.12. Exact pins in `requirements.txt`:

| Package | Version | Used for |
|---|---|---|
| numpy | 2.4.6 | float ops in the PU-bagging path |
| pandas | 3.0.5 | every table read and write |
| scipy | 1.18.0 | Fisher exact, Wilcoxon, permutation, FDR |
| scikit-learn | 1.9.0 | CV splitting and AUC |
| xgboost | 3.3.0 | base learner (`pu_target_model.py`) |
| pyarrow | 25.0.0 | parquet engine for `pu_model_matrix.parquet` |
| matplotlib | 3.11.0 | diagnostic plots only |

- `xgboost==3.3.0` fixes gene rank order and publishes wheels only for Python 3.12 and later.
- `src/analysis/common/_version_guard.py` runs `check_pins()` as step 1 of `make tables` and again when `pu_target_model.py` is imported, and aborts on a mismatch. `IGNITE_ALLOW_VERSION_MISMATCH=1` downgrades it to a warning; outputs from such a run should not be committed.
- `audit_coverage()` in the same module checks that every script constructing an `XGBClassifier` calls `check_pins()`.
- Unpinned, imported only by opt-in targets: `h5py`, `requests`, `openpyxl`. `pytest` backs `make test` only.

## R — stage 2, `make figures`

Requires R 4.4.1, parsed by `tools/check_r_versions.R` from `src/figures/R-requirements.txt`. Exact pins:

| Package | Version | Used for |
|---|---|---|
| ggplot2 | 4.0.3 | all figures |
| patchwork | 1.3.2 | multi-panel composition |
| dplyr | 1.2.1 | grouping in the genetic-data and specificity renderers |
| jsonlite | 2.0.0 | reads the `*_counts.json` / `*_stats.json` / `*_meta.json` tables |
| ragg | 1.5.2 | AGG PNG device |
| scales | 1.4.0 | axis and label formatting |
| ggrepel | 0.9.6 | label layout in `fig_discordance.R` |
| systemfonts | 1.2.1 | font matching for ragg |
| textshaping | 1.0.0 | glyph shaping for ragg |
| gridtext | 0.1.6 | layout engine behind ggtext |
| ggtext | 0.2.0 | `element_markdown()` / `geom_richtext()` in `fig_stat4_vignette.R` panel c |

Install exact versions:

```r
remotes::install_version("ggplot2", "4.0.3")
remotes::install_version("patchwork", "1.3.2")
remotes::install_version("dplyr", "1.2.1")
remotes::install_version("jsonlite", "2.0.0")
remotes::install_version("ragg", "1.5.2")
remotes::install_version("scales", "1.4.0")
remotes::install_version("ggrepel", "0.9.6")
remotes::install_version("systemfonts", "1.2.1")
remotes::install_version("textshaping", "1.0.0")
remotes::install_version("gridtext", "0.1.6")
remotes::install_version("ggtext", "0.2.0")
```

- `tools/setup_env.sh` installs in three passes; pass 3 fetches from the CRAN Archive via `remotes::install_version` to reach pins older than the current CRAN release.
- Native prerequisites, checked before any install: freetype (`ft2build.h`) for `ragg`, harfbuzz (`hb-ft.h`) for `textshaping`, and a working C++ toolchain. `IGNITE_CXX_SDK_FALLBACK=1` builds against the SDK's libc++.
- On macOS, `ggplot2`, `dplyr` and `ragg` need `type = "source"` to reach the pins; CRAN binaries are one patch behind. `gridtext` needs the CRAN `jpeg` package.
- System glyph libraries cannot be pinned by `src/figures/R-requirements.txt`. The verified build linked freetype 2.14.3, libpng 1.6.58, libtiff 4.7.2, jpeg-turbo 3.2.0 on macOS 26.5.1 arm64.

## Overrides

```
make tables  PY=/path/to/python3
make figures RSCRIPT=/path/to/Rscript
```

Both use `?=`, so an environment variable of the same name takes precedence. `PY` must point to the pinned environment.

# The figure_data contract

- `make tables` writes `figure_data/`; `make figures` reads only from it.
- `figure_data/` is committed, so stage 2 is standalone and needs no model fitting.
- A `make tables` run leaves the tree dirty under `figure_data/`. `make verify-tables` is the check that it reproduced.
- Numbered supplementary tables are packaged straight from `figure_data/` by `tools/package_supplementary_tables.py` and never enter `final_plots/`. Names and numbers live in `tools/supplementary_tables.tsv`.

# Stage-1 step order

`make tables` runs 30 Python steps. The order is load-bearing and the stage does not support `-j`. Paths below are relative to `src/analysis/`.

1. `common/_version_guard.py`
2. `evidence/make_enrichment_table.py`
3. `model/pu_target_model.py`
4. `model/attribution_decomposition.py`
5. `model/reshape_attribution_to_panels.py`
6. `model/compute_nb_significance.py`
7. `model/recompute_block_signal.py`
8. `model/make_panelD_recovery.py`
9. `model/build_score_by_group.py`
10. `evidence/build_genetic_data_tables.py`
11. `evidence/build_perturbational_fg_tables.py`
12. `vignette/build_vignette_panelB.py`
13. `vignette/build_vignette_panelD_knn.py`
14. `vignette/build_vignette_meta.py`
15. `specificity/build_specificity_matrix.py`
16. `specificity/compute_delong_specificity.py`
17. `specificity/build_leakage_controlled.py`
18. `validation/build_temporal_holdout.py`
19. `validation/build_heldout_trial_auc.py`
20. `tables/build_ranked_atlas_table.py`
21. `validation/build_immune_stages_data.py`
22. `tables/build_feature_dictionary.py`
23. `tables/build_drug_status_table.py`
24. `evidence/build_univariate_by_group.py`
25. `evidence/build_orthogonality_tables.py`
26. `model/build_coalition_tables.py`
27. `tables/build_supp_stat_sheets.py`
28. `vignette/build_knn_nearest_target.py`
29. `vignette/add_knn_approved_drugs.py`
30. `discordance/build_discordance.py`

Dependencies that fix the order:

| Step | Must follow | Because |
|---|---|---|
| `compute_nb_significance.py` | `reshape_attribution_to_panels.py` | rewrites `panelA_significance.csv` and `panelC_marginal.csv` in place |
| `build_vignette_meta.py` | `pu_target_model.py`, `build_vignette_panelD_knn.py` | reads `full_model_pu_scores.csv` and `vignette_knn_stats.json` |
| `build_ranked_atlas_table.py` | `pu_target_model.py` | joins `full_model_pu_scores.csv` onto the feature matrix |
| `build_immune_stages_data.py` | `build_ranked_atlas_table.py`, `make_panelD_recovery.py` | reads `ranked_atlas.csv` and `genetics_only_pu_scores.csv` |
| `build_feature_dictionary.py` | `pu_target_model.py`, `build_ranked_atlas_table.py` | reads `pu_feature_importance_stability.csv` and `ranked_atlas.csv` |
| `build_univariate_by_group.py`, `build_orthogonality_tables.py` | `build_ranked_atlas_table.py` | read `ranked_atlas.csv` |
| `build_leakage_controlled.py` | `pu_target_model.py` | reads `full_model_pu_scores.csv`; imports `delong_test` from `compute_delong_specificity.py` |
| `build_supp_stat_sheets.py` | steps 9, 17, 19, 21 and 13 | reshapes their JSON outputs |
| `add_knn_approved_drugs.py` | `build_knn_nearest_target.py` | annotates its `outputs/knn/` output |
| `build_discordance.py` | all of the above | runs last |

Stage-1 scripts write to two places: `figure_data/` (committed, checked by `verify-tables`) and `outputs/` (scratch, not committed, removed by `make clean`). `make tables` starts with `rm -rf logs` and writes per-step logs to `logs/NN_<script>.log`, with a summary in `logs/_summary.tsv`.

# Model configuration

- Base learner, seed and bagging config are defined once in `src/analysis/model/pu_target_model.py`: `SEED`, `T_BAG`, `base_learner()`, `pu_bag()`, `impute()`.
- Production ranking: `SEED = 4`, `T_BAG = 200`. Stage 2b re-scores every gene with 5 folds over the positives so no gene is scored in-sample.
- Coalition cross-validation (`attribution_decomposition.py`): all 8 coalitions of the three feature groups, 5-fold, averaged over 5 seeds.
- `recompute_block_signal.py` runs at `T_BAG = 60`.
- Feature groups, defined in `attribution_decomposition.py`: `GENETIC` (5) + `PERTURBATIONAL` (15) + `OBSERV` (4) = 24, the matrix feature count.
- Multiple-testing correction is applied within each pre-specified family: `holm()` for small pre-specified comparison sets, `bh()` for screens. Both live in `src/analysis/common/_stats.py`.

# Verification

| Target | Checks | Behaviour |
|---|---|---|
| `make check-versions` | Python and R stacks against the pins; `XGBClassifier` guard coverage | exit 1 on mismatch |
| `make verify-tables` | `figure_data/` byte-identical to the committed state | hard failure |
| `make verify-figures` | `final_plots/` against `tools/final_plots.sha256` | PNG/PDF differences reported; fatal for absent baseline files and for files the baseline does not list |

- `verify-tables` is the check for the numbers and is portable across machines. `verify-figures` reports PNG and PDF differences rather than gating on them, because `ragg` rasterizes through freetype, libpng and the installed fonts, which no lockfile pins, and every R PDF device stamps a wall-clock creation date.
- `verify-figures` sorts each baseline file into one of five classes before reporting: MISSING (fatal by default; `VERIFY_ALLOW_MISSING` downgrades it, set only by the `supplementary` target), RENDERED (gitignored `.png`/`.pdf`; soft), DOCS (`.md`; reported only), PUBLISHED (gitignored non-image build output; fatal), COMMITTED (tracked, externally sourced; fatal).
- 11 of the 98 tables in `figure_data/` belong to opt-in targets that `make tables` does not run — the 6 `discordance_*` caches, the 2 `label_permutation_*` files plus `scrambled_feature_control.json`, and the 2 `vignette_panelC_gwas_*` files. `verify-tables` confirms these are unchanged rather than re-derived.
- `tools/write_baseline.sh <py> <rscript>` regenerates `tools/final_plots.sha256` from a clean clone built with the pinned interpreters. Run it through the Makefile so R resolves packages from `.rlib/`.
- `tools/compare_submission.py` compares two assembled `final_outputs/` packages, normalising the wall-clock fields in `.xlsx` (`docProps/core.xml`) and `.pdf` (`/CreationDate`, `/ModDate`, `/ID`) before hashing.

# Opt-in targets

Six targets sit outside `make tables`, which stays offline and completes in minutes:

| Target | Time | Needs | Runs | Rebuilds |
|---|---|---|---|---|
| `make permutation-null` | ~8.5 h at `PERM_JOBS=11` | compute | `permutation/run_label_permutation.py`, `permutation/scrambled_feature_control.py` | `label_permutation_null.csv`, `label_permutation_pvalues.csv`, `scrambled_feature_control.json` |
| `make permutation-null-finalize` | seconds | — | `permutation/finalize_label_permutation.py` | p-values only, from the existing null |
| `make discordance-network` | ~10 min first run, ~1 min after | network, ~135 MB download | `discordance/fetch_discordance_network.py` | the 6 `discordance_*` cache tables |
| `make vignette-panelc` | seconds | network | `vignette/fetch_vignette_panelC_gwas.py` | `vignette_panelC_gwas_diseases.csv`, `vignette_panelC_gwas_meta.json` |
| `make gps-inputs` | seconds | network | `specificity/build_gps_per_gene.py` | `data/specificity/gps_per_gene.csv`, `gps_drug_labeled_genes.csv` |
| `make knn-signatures` | ~3-4 min | network | `vignette/build_knn_signatures.py`; then rerun `vignette/build_vignette_panelD_knn.py` as it prints | `data/perturbseq/knn/signatures_*.npy` |

- `permutation-null` checkpoints every 25 permutations into gitignored `outputs/permutation/`, so an interrupted run resumes. Smoke test with `make permutation-null PERM_N=8 PERM_SCR=2`.
- These targets write directly into `figure_data/` and have no dry-run mode. Revert a smoke run with `git checkout -- figure_data`.
- `make vignette-panelc` asserts Open Targets release 26.06 and aborts under any other release.
- `make discordance-network` downloads land in gitignored `data/raw/discordance/` and are reused on later runs; `make clean` does not remove them.
- `build_vignette_panelD_knn.py`, `build_knn_nearest_target.py` and `add_knn_approved_drugs.py` skip with a message when the signature matrices are absent, so `make tables` and `make verify-tables` complete on a clean clone. `src/analysis/common/_skip_ledger.py` prints which steps skipped at the end of the run.

# Platform

- Reference platform: Apple Silicon macOS (arm64). Verified under Python 3.12.13, R 4.4.1, macOS 26.5.1. There `make tables` reproduces all 98 tables byte-identically and `verify-tables` fails hard on any difference.
- Off the reference platform `verify-tables` reports which tables differ and exits successfully, and checks that tie a table to the committed run warn rather than stop the build. `src/analysis/common/_platform.py` holds the platform test they share.
- `vector_output.R` writes PDFs through the quartz device (macOS only). Elsewhere each renderer writes its PNG, logs that the PDF was skipped, and continues; a skipped render also deletes any older PDF of that figure. `IGNITE_SKIP_PDF=1` forces the same path on macOS.
- PNG bytes depend on the font and rasterization stack, which no lockfile pins, so they differ across platforms even when every pinned package matches. `verify-figures` treats this as a soft difference.
- The generated `final_outputs/README.md` states the platform it was built on.

# Running one step

The unit of work is one script. To rerun a single step as the build does:

```
PY=.venv/bin/python RSCRIPT=Rscript tools/run_step.sh py src/analysis/discordance/build_discordance.py
PY=.venv/bin/python RSCRIPT=Rscript tools/run_step.sh R  src/figures/fig_discordance.R
```

Or directly:

```
PYTHONPATH=src .venv/bin/python src/analysis/specificity/build_specificity_matrix.py
R_LIBS_USER=$PWD/.rlib Rscript src/figures/fig_specificity.R
```

# Conventions

- Python scripts locate the repository root by walking up from `__file__`, never from the current working directory.
- Modules import each other by bare name; a script needing a module from another folder puts that folder on `sys.path`. Run every script by path from the repo root, as the Makefile does.
- Counts are computed rather than hardcoded where possible.
- Build outputs are not committed; a clean clone regenerates them.
