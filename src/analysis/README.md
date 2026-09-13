# Analysis Layout and Files

## Layout

| Folder | Description |
|---|---|
| `common/` | shared helpers: `_version_guard`, `_stats`, `_platform`, `_skip_ledger`, `_gene_names` |
| `model/` | the PU model and its decomposition (Fig 4, Supplementary Table 6) |
| `evidence/` | genetic and perturbational evidence tables (Fig 2, Fig 3) |
| `validation/` | trial and temporal validation (Fig 5, Fig 6, Supplementary Table 7) |
| `specificity/` | specificity and leakage control (Fig 7, S3) |
| `vignette/` | STAT4 vignette and kNN signatures (Fig 8, Supplementary Tables 14-15) |
| `permutation/` | label-permutation and scrambled-feature controls (S1) |
| `discordance/` | discordance analysis (S2) |
| `tables/` | cross-cutting supplementary tables |

Modules import each other by bare name. A script that needs a module from another folder puts that folder on `sys.path` next to its own, so run every script by path from the repo root, as the Makefile does.

To re-derive the counts below, see [Re-deriving the counts](#re-deriving-the-counts).

---

## Files

### The reference implementation

`pu_target_model.py` defines `SEED`, `T_BAG`, `base_learner()`, `pu_bag()` and `impute()`, and runs the version guard at import time. Importing it inherits the config; a script that re-declares the learner instead can train under any installed xgboost. `_version_guard.audit_coverage()` requires every file constructing an `XGBClassifier` to call `check_pins()` or import a module that does.

### Inputs

| Path | Description |
|---|---|
| `data/perturbseq/pu/pu_model_matrix.parquet` | 19,502 × 28 = **24 features + 4 meta**. `pu_role`: 340 `P`, 727 `trial_heldout`, 18,435 `unlabeled` |
| `data/perturbseq/pu/pu_labels.csv` | the label table; the `pu_role` partition `build_perturbational_fg_tables.py` reads |

`scripts/build_pu_matrix.py`, which builds the parquet from the by-gene drug-target file, is not in this repository. The parquet is the committed starting point, so the feature set can only be consumed here, not re-derived.

### Scripts that train a PU model (10)

Grouped below by how each script obtains the config:

*Inherit from `pu_target_model` — config cannot drift:*

| Script | Description |
|---|---|
| `build_discordance.py` | `import pu_target_model as pum`; refits both arms at `pum.T_BAG` / `pum.SEED` |
| `build_heldout_trial_auc.py` | `import pu_target_model as pum`; re-scores the genetics-only comparator |
| `make_panelD_recovery.py` | loads `pum` via `importlib.spec_from_file_location` |

*Inherit from `attribution_decomposition` (which has its own learner — see below):*

| Script | Description |
|---|---|
| `compute_nb_significance.py` | calls `ad.cv_auc` over 4 models × 5 seeds × 5 folds, so the NB correction is scored identically to the panel it corrects |
| `scrambled_feature_control.py` | rides on `ad`'s learner and `ad.N_SEED` |
| `run_label_permutation.py` | uses `ad`'s groups and `cv_auc`, but **monkeypatches `ad.base_learner`** with its own `n_jobs=1` copy for thread-safety |

*Four of the five copies declare their own learner; `run_label_permutation` above is the fifth.*

| Script | Description |
|---|---|
| `pu_target_model.py` | the reference |
| `attribution_decomposition.py` | own `base_learner`, own inline bagging loop |
| `recompute_block_signal.py` | own `base_learner`, own `pu_bag_score` |
| `retrain_at_T0.py` | own `base_learner`, own `pu_bag`; a library, called by `build_temporal_holdout.py` |

## Settings

### Feature groups

`GENETIC` (5) + `PERTURBATIONAL` (15) + `OBSERV` (4) = **24**, the matrix feature count; the two sets are equal, not merely equinumerous. The canonical definition is in `attribution_decomposition.py` (`GENETIC`/`PERTURBATIONAL`/`OBSERV`, just below `META`). `GENETIC` is re-declared in `build_heldout_trial_auc.py`, `build_discordance.py`, and `make_panelD_recovery.py`, and, as `GENETICS_FEATS` in a different order, in `retrain_at_T0.py`.

---

## Re-deriving the counts

```bash
# matrix shape and label balance
python -c "import pandas as pd; m=pd.read_parquet('data/perturbseq/pu/pu_model_matrix.parquet'); print(m.shape); print(m.pu_role.value_counts())"

# every learner and bagging-loop copy
grep -rn "def base_learner\|def _single_thread_learner\|XGBClassifier(" --include="*.py" src
grep -rn "def pu_bag\|pu_bag(" --include="*.py" src

# feature groups still partition the matrix
grep -n "^GENETIC\|^PERTURBATIONAL\|^OBSERV" src/analysis/model/attribution_decomposition.py
```
