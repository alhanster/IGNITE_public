# PU model — file inventory and settings census

What trains a PU model, what consumes its scores, and every setting either of them depends on.

The PU model has wide reach in the pipeline: ten scripts train one, seven more read its scores, and together they feed every display item except Fig 2. Settings live in five `base_learner` definitions and four bagging loops rather than in one place, so a config drift in any copy changes published numbers with no error or traceback. `_version_guard.py` closes the same failure mode for the xgboost pin.

## Layout

| Folder | Scripts for |
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

Counts below are for orientation only. See [Re-deriving the counts](#re-deriving-the-counts) for the derivation, and `figure_data/PROVENANCE.md` for the note on integers that are not recomputed.

---

## 1. Files

### The reference implementation

`pu_target_model.py` defines `SEED`, `T_BAG`, `base_learner()`, `pu_bag()` and `impute()`, and runs the version guard at import time. Importing it inherits the config; a script that re-declares the learner instead can train under any installed xgboost. `_version_guard.audit_coverage()` requires every file constructing an `XGBClassifier` to call `check_pins()` or import a module that does.

### Inputs — committed

| Path | What |
|---|---|
| `data/perturbseq/pu/pu_model_matrix.parquet` | 19,502 × 28 = **24 features + 4 meta**. `pu_role`: 340 `P`, 727 `trial_heldout`, 18,435 `unlabeled` |
| `data/perturbseq/pu/pu_labels.csv` | the label table; the `pu_role` partition `build_perturbational_fg_tables.py` reads |

`scripts/build_pu_matrix.py`, which builds the parquet from the by-gene drug-target file, is not in this repository. The parquet is the committed starting point, so the feature set can only be consumed here, not re-derived.

### Outputs — under `outputs/model/`

**Nothing under `outputs/` is committed**; the directory is stage-1 scratch, recreated by `make tables`. `full_model_pu_scores.csv` is the artifact everything downstream reads. `genetics_only_pu_scores.csv`, written by `make_panelD_recovery.py`, is the genetics-only comparison arm behind Fig 4 panel e, using the same PU-bagging config restricted to the four `GENETIC` features.

`genetics_only_pu_scores.csv` also feeds `build_immune_stages_data.py`, which ranks it the same way as the full model, giving Fig 5 panels b and c their genetics-only arm. `make_panelD_recovery.py` is therefore an ordering constraint on the Fig 5 step as well as on Fig 4 panel e. `pu_role` is the same label column in both files, so the non-approved set is the same 19,162 genes for both arms; only the order changes.

The scorer writes exactly one file into the committed stage-1 contract:
`figure_data/vignette_score_hist.csv`.

### Scripts that TRAIN a PU model (10)

Grouped below by how each script obtains the config:

*Inherit from `pu_target_model` — config cannot drift:*

| Script | How |
|---|---|
| `build_discordance.py` | `import pu_target_model as pum`; refits both arms at `pum.T_BAG` / `pum.SEED` |
| `build_heldout_trial_auc.py` | `import pu_target_model as pum`; re-scores the genetics-only comparator |
| `make_panelD_recovery.py` | loads `pum` via `importlib.spec_from_file_location` |

*Inherit from `attribution_decomposition` (which has its own learner — see below):*

| Script | How |
|---|---|
| `compute_nb_significance.py` | calls `ad.cv_auc` over 4 models × 5 seeds × 5 folds, so the NB correction is scored identically to the panel it corrects |
| `scrambled_feature_control.py` | rides on `ad`'s learner and `ad.N_SEED` |
| `run_label_permutation.py` | uses `ad`'s groups and `cv_auc`, but **monkeypatches `ad.base_learner`** with its own `n_jobs=1` copy for thread-safety |

*Four of the five copies declare their own learner; `run_label_permutation` above is the fifth.*

| Script | Note |
|---|---|
| `pu_target_model.py` | the reference |
| `attribution_decomposition.py` | own `base_learner`, own inline bagging loop |
| `recompute_block_signal.py` | own `base_learner`, own `pu_bag_score` |
| `retrain_at_T0.py` | own `base_learner`, own `pu_bag`; a library, called by `build_temporal_holdout.py` |

### Scripts that CONSUME `full_model_pu_scores.csv` without refitting (8)

`build_score_by_group.py`, `build_immune_stages_data.py`, `build_specificity_matrix.py`,
`compute_delong_specificity.py`, `build_leakage_controlled.py`, `build_vignette_meta.py`,
`build_vignette_panelD_knn.py`, `build_knn_nearest_target.py`. `build_immune_stages_data.py`
also reads `genetics_only_pu_scores.csv`, ranking both files through one function that differs
between the arms only in the score column.

`retrain_at_T0.py`, `_version_guard.py` and `src/figures/fig_stat4_vignette.R` mention
`full_model_pu_scores.csv` in comments but do not read it. Panel a reads
`figure_data/vignette_score_hist.csv`, a committed trim of that ranking; every stage-2 script
reads `figure_data/` only.

### Figures downstream

Fig 4 (model analysis), Fig 5 (trial validation), Fig 6 (temporal holdout), Fig 7
(specificity), Fig 8 (STAT4 vignette), and S1-S3. The held-out in-trial AUC
statistics ship as the ROC and discrimination sheets of Supplementary Table 7, sourced from
`figure_data/`.

Fig 2 is the only display item with no PU dependency. Fig 3 depends on the PU labels but not
on the fitted model: `build_perturbational_fg_tables.py` partitions on `pu_role` from `pu_labels.csv`,
so changing `SEED` or `T_BAG` cannot move it.

---

## 2. Settings

### Production config — `pu_target_model.py`

| Setting | Value | Note |
|---|---|---|
| `SEED` | `4` | representative seed: the draw closest to the 5-seed mean across panels a–e. Was `0` until 2026-07-23 |
| `T_BAG` | `200` | PU-bagging iterations for the full model |
| `N_CV` / `N_REPEAT` | `5` / `10` | the CV loop bags at `T=40`, **not** `T_BAG` |
| CV fold seed | `SEED + r` | per repeat |
| CV bag seed | `SEED + r * 97` | per repeat |
| per-bag learner seed | `seed + t` | |
| `n_estimators` | `120` | |
| `max_depth` | `2` | `MODEL_MAX_DEPTH` env override; a parsimony choice, no worse than depth 3 on CV AUC or held-out recovery |
| `learning_rate` | `0.1` | |
| `subsample` / `colsample_bytree` | `0.8` / `0.8` | |
| `eval_metric` / `n_jobs` / `verbosity` | `logloss` / `4` / `0` | |
| imputation | median, fit per bag on that bag's training rows only | leakage-safe |
| pseudo-negatives | `rng.choice(U_idx, size=|P|, replace=True)` | positives stay in-bag every round; their scores are optimistic by construction |

### Per-site divergence

Every site agrees with the production learner on `n_estimators`, `max_depth`, `learning_rate`,
`subsample`, `colsample_bytree` and `eval_metric`. They differ in:

| Script | `T_BAG` | seed(s) | `n_jobs` | per-bag learner seed |
|---|---|---|---|---|
| `pu_target_model` | 200 | 4 | 4 | `seed + t` |
| `retrain_at_T0` | 200 | 4 | 4 | `seed + t` |
| `attribution_decomposition` | 30 | 0–4 | 4 | **`t`** |
| `recompute_block_signal` | 60 | 4; `SEEDS` 0–4 | 4 | **`t`** |
| `run_label_permutation` | 30 (via `ad`) | `10_000 + i` | **1** | `t` |
| `scrambled_feature_control` | 30 (via `ad`) | `700 + i` | 4 | `t` |

The differing `T_BAG` values reflect cost trade-offs: 200 where a gene-level ranking is
published, 30-60 where an aggregate metric is averaged over five seeds. `n_jobs=1` in
`run_label_permutation` is required, since it parallelises across permutations, not within a fit.

`run_label_permutation`'s per-permutation seed `RandomState(10_000 + i)` makes `PERM_JOBS`
change wall-clock only, never the numbers.

### Seeds that are NOT model seeds

Additional seeds: `build_temporal_holdout.SEED = 1234` (permutation and bootstrap; the model is seeded separately by `retrain_at_T0`) and `build_immune_stages_data.SEED = 1234` (B = 10,000; N_TOP = 50; CTRL_RANK_MIN = 200; `DEPTHS_C = (25, 50, 100)` for Fig 5 panel c). Control draw size follows shortlist depth, not `N_CTRL`; `N_CTRL = 50` remains only as a reported JSON field, matching the depth-50 cell of Supplementary Table 10. `RandomState` is re-seeded per (arm, depth), so depth 50 matches the single-depth result and panels b and c never share a bar height. `build_specificity_matrix` uses a 2,000-draw bootstrap; `build_discordance.N_NULL = 20` is the feature-test null, and the network null (500) is in the fetch script.

Makefile knobs: `PERM_N = 1000`, `PERM_JOBS = 11`, `PERM_SCR = 99`, `DISC_NULL = 500`.

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

# who reads the scores (drop comment-only hits by eye)
grep -rn "full_model_pu_scores" --include="*.py" --include="*.R" src

# feature groups still partition the matrix
grep -n "^GENETIC\|^PERTURBATIONAL\|^OBSERV" src/analysis/model/attribution_decomposition.py
```
