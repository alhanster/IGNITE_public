# src/analysis — stage 1

30 Python steps: fit the models, run the tests, write `figure_data/`.

- Run them in the order the Makefile's `tables:` target gives. The order is load-bearing and the stage does not support `-j`.
- The step list and its dependency constraints are in `REPRODUCIBILITY.md`.
- Modules import each other by bare name. A script needing a module from another folder puts that folder on `sys.path` next to its own, so run every script by path from the repo root, as the Makefile does.

# Layout

| Folder | Contents |
|---|---|
| `common/` | shared helpers: `_version_guard`, `_stats`, `_platform`, `_skip_ledger`, `_gene_names` |
| `model/` | the PU model and its decomposition (Fig 4, Supp Table 6) |
| `evidence/` | genetic and perturbational evidence tables (Fig 2, Fig 3, Supp Tables 4-5) |
| `validation/` | trial and temporal validation (Fig 5, Fig 6, Supp Tables 7, 10, 12) |
| `specificity/` | specificity and leakage control (Fig 7, S3, Supp Tables 11, 13) |
| `vignette/` | STAT4 vignette and kNN signatures (Fig 8, Supp Tables 14-15) |
| `permutation/` | label-permutation and scrambled-feature controls (S1, Supp Table 8) |
| `discordance/` | discordance analysis (S2, Supp Table 9) |
| `tables/` | cross-cutting supplementary tables (Supp Tables 1-3, and the flattened stats sheets) |

# Inputs

| Path | Contents |
|---|---|
| `data/perturbseq/pu/pu_model_matrix.parquet` | 19,502 x 28 = 24 features + 4 meta. `pu_role`: 340 `P`, 727 `trial_heldout`, 18,435 `unlabeled` |
| `data/perturbseq/pu/pu_labels.csv` | the label table; the `pu_role` partition `build_perturbational_fg_tables.py` reads |

`scripts/build_pu_matrix.py`, which builds the parquet from the by-gene drug-target file, is not in this repository. The parquet is the committed starting point.

Re-derive the shape and label balance:

```bash
python -c "import pandas as pd; m=pd.read_parquet('data/perturbseq/pu/pu_model_matrix.parquet'); print(m.shape); print(m.pu_role.value_counts())"
```

# Outputs

- `figure_data/` — the committed contract read by stage 2 and checked by `make verify-tables`. One line per file in `figure_data/PROVENANCE.md`.
- `outputs/` — scratch that no figure reads. Not committed; removed by `make clean`.

# The PU model

- `model/pu_target_model.py` is the reference implementation: it defines `SEED`, `T_BAG`, `base_learner()`, `pu_bag()` and `impute()`, and runs the version guard at import.
- Importing it inherits the config; a script that re-declares the learner instead can train under any installed xgboost.
- `_version_guard.audit_coverage()` requires every file constructing an `XGBClassifier` to call `check_pins()` or import a module that does.

Ten scripts train a PU model, grouped by how each obtains its config.

Inherit from `pu_target_model`:

| Script | Notes |
|---|---|
| `discordance/build_discordance.py` | `import pu_target_model as pum`; refits both arms at `pum.T_BAG` / `pum.SEED` |
| `validation/build_heldout_trial_auc.py` | `import pu_target_model as pum`; re-scores the genetics-only comparator |
| `model/make_panelD_recovery.py` | loads `pum` via `importlib.spec_from_file_location` |

Inherit from `attribution_decomposition`:

| Script | Notes |
|---|---|
| `model/compute_nb_significance.py` | calls `ad.cv_auc` over 4 models x 5 seeds x 5 folds |
| `permutation/scrambled_feature_control.py` | rides on `ad`'s learner and `ad.N_SEED` |
| `permutation/run_label_permutation.py` | uses `ad`'s groups and `cv_auc`, but monkeypatches `ad.base_learner` with its own `n_jobs=1` copy for thread safety |

Declare their own learner:

| Script | Notes |
|---|---|
| `model/pu_target_model.py` | the reference |
| `model/attribution_decomposition.py` | own `base_learner`, own inline bagging loop |
| `model/recompute_block_signal.py` | own `base_learner`, own `pu_bag_score` |
| `validation/retrain_at_T0.py` | own `base_learner`, own `pu_bag`; a library, called by `build_temporal_holdout.py` |

# Feature groups

- `GENETIC` (5) + `PERTURBATIONAL` (15) + `OBSERV` (4) = 24, the matrix feature count.
- The canonical definition is in `model/attribution_decomposition.py`, below `META`.
- `GENETIC` is re-declared in `build_heldout_trial_auc.py`, `build_discordance.py` and `make_panelD_recovery.py`, and as `GENETICS_FEATS` in a different order in `retrain_at_T0.py`.

# Files `make tables` does not run

| File | What it is |
|---|---|
| `common/_stats.py` | `holm()` and `bh()`, the shared multiple-testing helpers |
| `common/_platform.py` | the reference-platform test shared by the committed-run checks |
| `common/_skip_ledger.py` | prints which steps skipped at the end of a run |
| `common/_gene_names.py` | gene-symbol helpers, including `load_hgnc` |
| `validation/retrain_at_T0.py` | library; the T0 retrain called by `build_temporal_holdout.py` |
| `validation/recovery_eval.py` | library; recovery@k with its bootstrap CI and permutation null |
| `permutation/_label_permutation.py` | library; shared paths and the p-value table |
| `permutation/run_label_permutation.py` | opt-in, `make permutation-null` |
| `permutation/scrambled_feature_control.py` | opt-in, `make permutation-null` |
| `permutation/finalize_label_permutation.py` | opt-in, `make permutation-null-finalize` |
| `discordance/fetch_discordance_network.py` | opt-in, `make discordance-network` |
| `vignette/fetch_vignette_panelC_gwas.py` | opt-in, `make vignette-panelc` |
| `vignette/build_knn_signatures.py` | opt-in, `make knn-signatures` |
| `specificity/build_gps_per_gene.py` | opt-in, `make gps-inputs` |

`common/_version_guard.py` is step 1 of `make tables` rather than a library.

Re-derive:

```bash
# every learner and bagging-loop copy
grep -rn "def base_learner\|def _single_thread_learner\|XGBClassifier(" --include="*.py" src
grep -rn "def pu_bag\|pu_bag(" --include="*.py" src

# feature groups still partition the matrix
grep -n "^GENETIC\|^PERTURBATIONAL\|^OBSERV" src/analysis/model/attribution_decomposition.py
```
