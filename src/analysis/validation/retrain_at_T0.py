#!/usr/bin/env python
"""Temporal-holdout retrain of the PU-bagging immune-target model at a past vintage T0.

Library entry point: run_T0() takes the label sets as arguments, defined by
build_temporal_holdout.build_split(). Trains a full model (24 features) and a
5-feature genetics-only comparator over the same 19,502 genes.

See REPRODUCIBILITY.md, Temporal holdout (T0) labels, for the split
design and the feature-timestamping caveat.
"""
import os, json
import numpy as np
import pandas as pd
from xgboost import XGBClassifier

# Aborts on an xgboost version mismatch, since a mismatch silently changes gene-level results.
import sys as _sys, os as _os
_r = _os.path.abspath(__file__)
while _r != _os.path.dirname(_r) and not _os.path.isfile(_os.path.join(_r, "Makefile")):
    _r = _os.path.dirname(_r)
_sys.path.insert(0, _os.path.join(_r, "src", "analysis", "common"))
from _version_guard import check_pins  # noqa: E402
check_pins()
# -----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MATRIX = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
META = ["gene", "pu_role", "furthest_stage", "crossdonor_correlation_mean"]
# Same genetic block as attribution_decomposition.py, build_discordance.py,
# build_heldout_trial_auc.py, and make_panelD_recovery.py, but in this script's own
# column order (gwas_score before IEI); colsample_bytree makes order load-bearing, so
# it is left as-is rather than aligned. See REPRODUCIBILITY.md, Cross-validation and seeds.
GENETICS_FEATS = ["lof.oe_ci.upper", "mis.z_score", "gwas_score", "IEI", "gene_burden_score"]
# SEED=4 fixes pu_target_model's medoid seed, keeping Fig 6 and gene-level
# claims aligned with outputs/model/full_model_pu_scores.csv. See
# REPRODUCIBILITY.md, Cross-validation and seeds.
SEED = 4
T_BAG = 200


def base_learner(seed):
    # Base_learner config must match pu_target_model.base_learner; MODEL_MAX_DEPTH overrides depth.
    return XGBClassifier(n_estimators=120,
                         max_depth=int(os.environ.get("MODEL_MAX_DEPTH", "2")),
                         learning_rate=0.1,
                         subsample=0.8, colsample_bytree=0.8, eval_metric="logloss",
                         n_jobs=4, verbosity=0, random_state=seed)


def impute(train_X, *apply_X):
    med = np.nanmedian(train_X, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    out = [np.where(np.isnan(train_X), med, train_X)]
    for A in apply_X:
        out.append(np.where(np.isnan(A), med, A))
    return out


def pu_bag(X_all, P_pos_idx, U_idx, T, seed):
    """PU-bagging; mean out-of-bag score per row (nan if the row was never held out across bags). Positive-labeled rows are kept in-bag."""
    rng = np.random.RandomState(seed)
    n = X_all.shape[0]
    ssum = np.zeros(n); scount = np.zeros(n)
    for t in range(T):
        neg = rng.choice(U_idx, size=len(P_pos_idx), replace=True)
        tr = np.concatenate([P_pos_idx, neg])
        y = np.r_[np.ones(len(P_pos_idx)), np.zeros(len(neg))]
        Xtr, Xall_i = impute(X_all[tr], X_all)
        clf = base_learner(seed + t)
        clf.fit(Xtr, y)
        p = clf.predict_proba(Xall_i)[:, 1]
        oob = np.ones(n, bool); oob[neg] = False
        ssum[oob] += p[oob]; scount[oob] += 1
    return np.where(scount > 0, ssum / np.maximum(scount, 1), np.nan)


def run_T0(pos_at_t0_genes, emergent_genes, T0, feature_subset=None, tag="full",
           holdout_genes=None):
    """Retrain at T0.

    pos_at_t0_genes, emergent_genes: sets of gene symbols.
    holdout_genes: genes known-in-trial-but-not-approved at T0; excluded from
    the pseudo-negative pool but not part of the emergent test set.

    Returns a per-gene DataFrame with score and rank_pctile among the scorable pool.
    """
    m = pd.read_parquet(MATRIX)
    all_feats = [c for c in m.columns if c not in META]
    feats = feature_subset if feature_subset else all_feats
    X = m[feats].apply(pd.to_numeric, errors="coerce").values.astype(float)
    gene = m["gene"].values

    holdout_genes = holdout_genes or set()
    is_P = np.isin(gene, list(pos_at_t0_genes))
    is_E = np.isin(gene, list(emergent_genes))
    is_H = np.isin(gene, list(holdout_genes))
    # Emergent and holdout genes are scored but never drawn as pseudo-negatives, matching how
    # trial_heldout is treated in the production model, which is what keeps the test set held out.
    is_U = ~is_P & ~is_E & ~is_H
    P_idx = np.where(is_P)[0]; U_idx = np.where(is_U)[0]

    score = pu_bag(X, P_idx, U_idx, T=T_BAG, seed=SEED)
    # Rank percentile is computed over the discovery pool only, excluding positives (in-bag,
    # would rank optimistically) and held-out-known genes (already in clinic).
    pool = is_E | is_U
    s = score.copy()
    rp = pd.Series(np.where(pool, s, np.nan)).rank(pct=True).values
    out = pd.DataFrame({"gene": gene,
                        f"score_{tag}": score,
                        f"rank_pctile_{tag}": rp,
                        "is_pos_at_T0": is_P, "is_emergent": is_E,
                        "is_holdout": is_H, "in_pool": pool})
    return out, dict(n_P=int(is_P.sum()), n_E=int(is_E.sum()), n_H=int(is_H.sum()),
                     n_U=int(is_U.sum()), n_feat=len(feats))


if __name__ == "__main__":
    print("module ready; call run_T0() from driver")
