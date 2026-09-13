#!/usr/bin/env python3
"""Figure 4 panel f: perturb-seq block signal by genetic-prior stratum.
See REPRODUCIBILITY.md for the feature groups.
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier

# Aborts on an unpinned xgboost version: a mismatch can silently change gene-level results (see _version_guard.py).
import sys as _sys, os as _os
_r = _os.path.abspath(__file__)
while _r != _os.path.dirname(_r) and not _os.path.isfile(_os.path.join(_r, "Makefile")):
    _r = _os.path.dirname(_r)
_sys.path.insert(0, _os.path.join(_r, "src", "analysis", "common"))
from _version_guard import check_pins  # noqa: E402
check_pins()
# -----------------------------------------------------------------------------

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
PU_DIR = os.path.join(ROOT, "data", "perturbseq", "pu")
OUTD   = os.path.join(ROOT, "figure_data")
os.makedirs(OUTD, exist_ok=True)
MATRIX = os.path.join(PU_DIR, "pu_model_matrix.parquet")
OUT    = os.path.join(OUTD, "panelE_block_signal.csv")

# Panel f aggregation and bag count: see REPRODUCIBILITY.md.
SEED, T_BAG = 4, 60
SEEDS = [0, 1, 2, 3, 4]
PRIOR = ["lof.oe_ci.upper", "mis.z_score", "IEI", "gwas_score", "gene_burden_score"]
META  = ["gene", "pu_role", "furthest_stage", "crossdonor_correlation_mean"]

# Includes only perturbational blocks; coverage flags are listed standalone rather than nested under their block.
BLOCKS = {
    "Trans-regulatory\nresidual": ["expected_n_regulators_residuals"],
    "Regulatory\nburden":         ["reg_burden_sig_Rest", "reg_burden_sig_Stim8hr",
                                    "reg_burden_sig_Stim48hr"],
    "Cytokine\nregulation":       ["n_sig_regulated_cytokines_Rest",
                                    "n_sig_regulated_cytokines_Stim8hr",
                                    "n_sig_regulated_cytokines_Stim48hr"],
    "Cytokine-receptor\nregulation": ["n_sig_regulated_cytokine_receptors_Rest",
                                       "n_sig_regulated_cytokine_receptors_Stim8hr",
                                       "n_sig_regulated_cytokine_receptors_Stim48hr"],
    "Polarization\nknockdown-coefficients": ["polar_coef_rank_Rest", "polar_coef_rank_Stim8hr",
                                              "polar_coef_rank_Stim48hr", "polar_rank_range"],
    # Feature matrix carries one combined coverage flag, not separate has_residual_burden, has_polarization, and known_polar_regulator flags.
    "Coverage\nflags":            ["has_cytokine"],
}


_MAX_DEPTH = int(os.environ.get("MODEL_MAX_DEPTH", "2"))   # final config = 2 (2026-07-23)


def base_learner(seed):
    # max_depth must match pu_target_model.base_learner, or this panel describes a different model than the published ranking. MODEL_MAX_DEPTH overrides only for depth-diff runs.
    return XGBClassifier(n_estimators=120, max_depth=_MAX_DEPTH, learning_rate=0.1, subsample=0.8,
                         colsample_bytree=0.8, eval_metric="logloss", n_jobs=4,
                         verbosity=0, random_state=seed)


def pu_bag_score(m, cols, Pidx, Uidx, T=T_BAG, seed=0):
    """Returns the out-of-bag PU-bagging score for each row; NaN where a row was never out-of-bag. Missing values are imputed with a leakage-safe median."""
    X = m[cols].apply(pd.to_numeric, errors="coerce").values.astype(float)
    rng = np.random.RandomState(seed)
    ssum = np.zeros(len(m)); scnt = np.zeros(len(m))
    for t in range(T):
        neg = rng.choice(Uidx, len(Pidx), replace=True)
        tri = np.concatenate([Pidx, neg])
        yt = np.r_[np.ones(len(Pidx)), np.zeros(len(neg))]
        med = np.nanmedian(X[tri], axis=0); med = np.where(np.isnan(med), 0.0, med)
        Xtr = np.where(np.isnan(X[tri]), med, X[tri])
        Xall = np.where(np.isnan(X), med, X)
        c = base_learner(t).fit(Xtr, yt)
        p = c.predict_proba(Xall)[:, 1]
        oob = np.ones(len(m), bool); oob[neg] = False
        ssum[oob] += p[oob]; scnt[oob] += 1
    return np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan)


def stratum_auc(score, mask, is_T, is_U):
    t = score[is_T & mask]; u = score[is_U & mask]
    y = np.r_[np.ones(len(t)), np.zeros(len(u))]; s = np.r_[t, u]
    v = ~np.isnan(s)
    return roc_auc_score(y[v], s[v])


def main():
    m = pd.read_parquet(MATRIX)
    # Sanity check: every block column is present.
    missing = [c for cols in BLOCKS.values() for c in cols if c not in m.columns]
    assert not missing, f"missing columns: {missing}"

    role = m.pu_role.values
    is_U, is_T = role == "unlabeled", role == "trial_heldout"
    Pidx, Uidx = np.where(role == "P")[0], np.where(is_U)[0]

    # Prior-only strata formed by median split of the unlabeled prior score.
    prior_score = pu_bag_score(m, PRIOR, Pidx, Uidx, T=100, seed=SEED)
    med = np.nanmedian(prior_score[is_U])
    weak = prior_score < med
    strong = ~weak

    # Prior-only baseline AUC per stratum computed as a 5-seed average.
    scp_avg = np.nanmean([pu_bag_score(m, PRIOR, Pidx, Uidx, seed=s)
                          for s in SEEDS], axis=0)
    base_w = stratum_auc(scp_avg, weak, is_T, is_U)
    base_s = stratum_auc(scp_avg, strong, is_T, is_U)

    rows = []
    for bname, bcols in BLOCKS.items():
        sc = np.nanmean([pu_bag_score(m, PRIOR + bcols, Pidx, Uidx, seed=s)
                         for s in SEEDS], axis=0)
        wl = stratum_auc(sc, weak, is_T, is_U) - base_w
        sl = stratum_auc(sc, strong, is_T, is_U) - base_s
        rows.append({"block": bname.replace("\n", " "), "weak_lift": wl, "strong_lift": sl})
        print(f"{bname.replace(chr(10),' '):32s} weak={wl:+.4f}  strong={sl:+.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nbaseline AUC: weak={base_w:.4f} strong={base_s:.4f}")
    print(f"n: trial_weak={int((is_T&weak).sum())} trial_strong={int((is_T&strong).sum())} "
          f"unlab_weak={int((is_U&weak).sum())} unlab_strong={int((is_U&strong).sum())}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
