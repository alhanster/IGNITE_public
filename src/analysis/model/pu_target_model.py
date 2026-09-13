#!/usr/bin/env python
"""PU-learning model that ranks immune drug-target candidates from perturb-seq and human-genetics features.

Labels follow a positive-unlabeled framing, since non-targets may be undiscovered targets rather than confirmed negatives. P consists of 340 genes with furthest_stage equal to Approved. trial_heldout consists of 727 Phase 1-3 genes, never trained on, used only to validate recovery and the Phase 3 > 2 > 1 dose-response gradient. The unlabeled pool (U) is the remaining 18,435 genes, the source of pseudo-negatives and the discovery space.

The model applies PU-bagging (Mordelet and Vert, 2014): each iteration draws a bootstrap sample of U as pseudo-negatives, trains an XGBoost classifier on P versus the pseudo-negatives, and scores the out-of-bag genes; a gene's final score is the mean over the iterations in which it was out-of-bag.

Metrics: recall@k of held-out trial targets is the headline result. The Phase 3 > 2 > 1 median-rank gradient validates dose-response. The Lee-Liu PU score guides tuning. ROC-AUC, with unlabeled treated as negative, serves model comparison only, since true positives hidden in U deflate its absolute value. The Elkan-Noto class prior estimates the number of undiscovered targets in U.

crossdonor_correlation_mean is joined to the output ranking as a post-hoc donor-reproducibility column and is not a model feature.

Requires scripts/build_pu_matrix.py to run first, producing pu_labels.csv and pu_model_matrix.parquet from the by-gene drug-target file.
"""
import os
import sys
import numpy as np
import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import roc_auc_score

# See REPRODUCIBILITY.md.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from _version_guard import check_pins  # noqa: E402
check_pins()
# -----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))   # evaluation/
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"   # repo root
EVAL = os.path.join(ROOT, "outputs", "model")   # scores/reports
PU_DIR = os.path.join(ROOT, "data", "perturbseq", "pu")   # input matrix
os.makedirs(EVAL, exist_ok=True)
MATRIX = os.path.join(PU_DIR, "pu_model_matrix.parquet")
OUT_RANK = os.path.join(EVAL, "full_model_pu_scores.csv")
OUT_REPORT = os.path.join(EVAL, "pu_model_report.txt")
OUT_IMP = os.path.join(ROOT, "outputs", "model", "pu_feature_importance.png")
FIGDATA = os.path.join(ROOT, "figure_data")
OUT_SCORE_HIST = os.path.join(FIGDATA, "vignette_score_hist.csv")

SEED = 4          # 2026-07-23: representative seed (medoid: draw closest to the 5-seed
                  # Mean across panels a-e (see seed_consistency_map).
T_BAG = 200            # PU-bagging iterations
N_CV = 5               # stratified folds over positives
N_REPEAT = 10          # CV repeats
META = ["gene", "pu_role", "furthest_stage", "crossdonor_correlation_mean"]


def base_learner(seed):
    # See REPRODUCIBILITY.md.
    return XGBClassifier(n_estimators=120, max_depth=int(os.environ.get("MODEL_MAX_DEPTH", "2")),
                         learning_rate=0.1,
                         subsample=0.8, colsample_bytree=0.8, eval_metric="logloss",
                         n_jobs=4, verbosity=0, random_state=seed)


def impute(train_X, *apply_X):
    """Median imputation; medians are fit on train_X only, to avoid leakage."""
    med = np.nanmedian(train_X, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    out = [np.where(np.isnan(train_X), med, train_X)]
    for A in apply_X:
        out.append(np.where(np.isnan(A), med, A))
    return out


def pu_bag(X_all, P_pos_idx, U_idx, feature_cols, T, seed, score_rows=None):
    """Runs PU-bagging. Returns the mean out-of-bag score for each row (NaN if never out-of-bag) and accumulated feature importances. Imputation medians are refit on each iteration's training rows only. score_rows selects which rows to score; defaults to all."""
    rng = np.random.RandomState(seed)
    n = X_all.shape[0]
    if score_rows is None:
        score_rows = np.arange(n)
    ssum = np.zeros(n); scount = np.zeros(n)
    imp_mat = np.zeros((T, len(feature_cols)))
    for t in range(T):
        neg = rng.choice(U_idx, size=len(P_pos_idx), replace=True)
        tr = np.concatenate([P_pos_idx, neg])
        y = np.r_[np.ones(len(P_pos_idx)), np.zeros(len(neg))]
        Xtr, Xall_i = impute(X_all[tr], X_all)
        clf = base_learner(seed + t)
        clf.fit(Xtr, y)
        p = clf.predict_proba(Xall_i)[:, 1]
        imp_mat[t] = clf.feature_importances_
        oob = np.ones(n, bool)
        oob[neg] = False           # rows used as pseudo-neg this round are not OOB
        # In-bag positive scores here are optimistic; Elkan-Noto prior uses CV held-out scores instead.
        ssum[oob] += p[oob]; scount[oob] += 1
    score = np.where(scount > 0, ssum / np.maximum(scount, 1), np.nan)
    return score, imp_mat.mean(0), imp_mat.std(0)


def recall_at_k(score, target_mask, unlabeled_mask, ks=(50, 100, 200, 500)):
    """Among the top-k genes ranked from target and unlabeled genes, recall_at_k reports the count or fraction that are true targets."""
    pool = target_mask | unlabeled_mask
    s = score.copy(); s[~pool] = -np.inf
    order = np.argsort(-s)
    out = {}
    n_tgt = int(target_mask.sum())
    for k in ks:
        topk = order[:k]
        hit = int(target_mask[topk].sum())
        out[k] = {"hits": hit, "precision": hit / k, "recall": hit / n_tgt}
    return out


def main():
    m = pd.read_parquet(MATRIX)
    feature_cols = [c for c in m.columns if c not in META]
    X = m[feature_cols].apply(pd.to_numeric, errors="coerce").values.astype(float)
    role = m["pu_role"].values
    is_P = role == "P"; is_U = role == "unlabeled"; is_T = role == "trial_heldout"
    P_idx = np.where(is_P)[0]; U_idx = np.where(is_U)[0]
    rep = []

    # Step 1: CV over positives measures recovery; also evaluates held-out trial genes at each config.
    skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=SEED)
    aucs, ll_scores, heldout_P = [], [], []
    for r in range(N_REPEAT):
        skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=SEED + r)
        for tr_p, te_p in skf.split(P_idx, np.ones(len(P_idx))):
            train_P = P_idx[tr_p]; test_P = P_idx[te_p]
            score, _, _ = pu_bag(X, train_P, U_idx, feature_cols, T=40, seed=SEED + r * 97)
            # ROC-AUC (held-out positives vs unlabeled) is for relative model comparison only, not absolute performance.
            mask = np.zeros(len(m), bool); mask[test_P] = True; mask[is_U] = True
            y = np.zeros(len(m)); y[test_P] = 1
            valid = mask & ~np.isnan(score)
            aucs.append(roc_auc_score(y[valid], score[valid]))
            heldout_P.append(np.nanmean(score[test_P]))
            # Lee-Liu score: recall of held-out positives at threshold = predicted-positive fraction of pool.
            thr = np.nanmedian(score[test_P])   # operating point
            pool = valid
            recall = np.mean(score[test_P] >= thr)
            ppp = np.mean(score[pool] >= thr)
            ll_scores.append(recall ** 2 / ppp if ppp > 0 else 0.0)
    rep.append(f"CV ROC-AUC (held-out approved vs unlabeled): "
               f"{np.mean(aucs):.3f} +/- {np.std(aucs):.3f}  (n={len(aucs)})")
    rep.append(f"Lee-Liu PU score (recall^2 / P(pred pos)): "
               f"{np.mean(ll_scores):.3f} +/- {np.std(ll_scores):.3f}")

    # ---- (2) Full model: trains on all approved genes and scores every gene ----
    score, imp, imp_sd = pu_bag(X, P_idx, U_idx, feature_cols, T=T_BAG, seed=SEED)
    rank_pct = pd.Series(score).rank(pct=True).values

    # ---- (2b) Cross-fitted scores: every gene scored by folds that did not train on it ----
    # Positives are in-bag on every iteration of the full fit, so their full-model scores are
    # in-sample. Five folds over P give each positive a score from the one fold where it was
    # held out; non-positives are credited in all five folds. Used for the drug-status group
    # comparison (Fig 4d, Supp Table 7), never for the published ranking.
    cv_sum = np.zeros(len(m)); cv_cnt = np.zeros(len(m))
    for tr_p, te_p in KFold(N_CV, shuffle=True, random_state=SEED).split(P_idx):
        s_cf, _, _ = pu_bag(X, P_idx[tr_p], U_idx, feature_cols, T=T_BAG, seed=SEED)
        credit = ~np.isnan(s_cf)
        credit[P_idx[tr_p]] = False          # fold's training positives are in-bag
        cv_sum[credit] += s_cf[credit]; cv_cnt[credit] += 1
    score_cv = np.where(cv_cnt > 0, cv_sum / np.maximum(cv_cnt, 1), np.nan)
    rank_pct_cv = pd.Series(score_cv).rank(pct=True).values
    assert np.isfinite(score_cv).all(), "cross-fit left a gene unscored"
    assert (cv_cnt[is_P] == 1).all() and (cv_cnt[~is_P] == N_CV).all(), \
        "cross-fit credit must be 1 fold per positive and N_CV folds per other gene"

    # ---- (3) Held-out trial recovery: trial genes never included in training, a genuine out-of-sample check ----
    rec = recall_at_k(score, is_T, is_U)
    rep.append("\nHeld-out trial-target recovery (top-k among trial ∪ unlabeled):")
    for k, v in rec.items():
        rep.append(f"  top-{k:4d}: {v['hits']:3d} trial hits  "
                   f"precision={v['precision']:.3f}  recall={v['recall']:.3f}")
    rep.append(f"  median rank-pctile: trial={np.median(rank_pct[is_T]):.3f}  "
               f"unlabeled={np.median(rank_pct[is_U]):.3f}  approved={np.median(rank_pct[is_P]):.3f}")
    rep.append(f"  median rank-pctile, cross-fitted: trial={np.median(rank_pct_cv[is_T]):.3f}  "
               f"unlabeled={np.median(rank_pct_cv[is_U]):.3f}  approved={np.median(rank_pct_cv[is_P]):.3f}")

    # ---- (4) Phase ladder dose-response ----
    rep.append("\nPhase-ladder gradient (median score-pctile by furthest_stage):")
    for stg in ["Approved", "Phase 3", "Phase 2", "Phase 1", "Unknown"]:
        sel = m["furthest_stage"].values == stg
        if sel.sum():
            rep.append(f"  {stg:9s} n={int(sel.sum()):4d}  median_pctile={np.median(rank_pct[sel]):.3f}")

    # ---- (5) Elkan-Noto class prior estimate; see REPRODUCIBILITY.md ----
    c_hat = float(np.nanmean(heldout_P))
    est_pos_frac = np.nanmean(score[is_U]) / c_hat if c_hat > 0 else np.nan
    est_pos_frac = float(np.clip(est_pos_frac, 0, 1))
    rep.append(f"\nElkan-Noto class prior (REPORTED WITH CAVEAT):")
    rep.append(f"  c_hat (mean score on CV held-out positives) = {c_hat:.3f}")
    rep.append(f"  implied positive fraction in U ~ {est_pos_frac:.3f} "
               f"(~{int(est_pos_frac * is_U.sum())} of {int(is_U.sum())})")
    rep.append(f"  *** NOT CREDIBLE as a literal count. c_hat ~ 0.5 means approved targets")
    rep.append(f"  are not cleanly separable from U at the gene-feature level, which violates")
    rep.append(f"  the Elkan-Noto well-separated-positive assumption. The prior is not")
    rep.append(f"  identifiable here; trust the recovery metrics, not this number. ***")

    # ---- (6) Output ranking: attaches a crossdonor confidence value to each ranked candidate ----
    out = m[["gene"]].copy()
    out["pu_score"] = score
    out["rank_pctile"] = rank_pct
    out["pu_score_cv"] = score_cv          # group comparison only; the ranking stays on pu_score
    out["rank_pctile_cv"] = rank_pct_cv
    out["pu_role"] = m["pu_role"].values
    out["furthest_stage"] = m["furthest_stage"].values
    out["crossdonor_confidence"] = m["crossdonor_correlation_mean"].values
    out = out.sort_values("pu_score", ascending=False, na_position="last").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out.to_csv(OUT_RANK, index=False)

    # See REPRODUCIBILITY.md.
    os.makedirs(FIGDATA, exist_ok=True)
    out[["gene", "pu_score"]].to_csv(OUT_SCORE_HIST, index=False)

    # Top-ranked unlabeled candidates (not approved positives), each with a donor-vetting flag.
    nov = out[out["pu_role"] == "unlabeled"].head(50).copy()
    cd = nov["crossdonor_confidence"]
    nov["donor_vetted"] = np.where(cd.notna(),
                                   np.where(cd >= 0.3, "yes", "weak"), "not_measured")
    nov.to_csv(os.path.join(EVAL, "top_novel_candidates.csv"), index=False)
    top100 = out.head(100)["pu_role"].value_counts().to_dict()
    rep.append(f"\nTop-100 composition: {top100}")
    rep.append("Top novel (unlabeled) candidates: "
               + ", ".join(nov.head(10)["gene"].tolist()))

    # feature importance & stability table
    imp_df = pd.DataFrame({"feature": feature_cols, "imp_mean": imp, "imp_sd": imp_sd})
    imp_df = imp_df.sort_values("imp_mean", ascending=False)
    imp_df.to_csv(os.path.join(EVAL, "pu_feature_importance_stability.csv"), index=False)

    # ---- (7) Feature importance & stability figure (mean +/- SD across bags) ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    top = imp_df.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    ax.barh(top["feature"], top["imp_mean"], xerr=top["imp_sd"], color="#56B4E9",
            edgecolor="black", linewidth=0.4,
            error_kw=dict(ecolor="#444444", lw=0.8, capsize=2))
    ax.set_xlabel("Mean XGBoost feature importance +/- SD across PU-bagging iterations (T=%d)" % T_BAG)
    ax.set_title("PU model feature importance & stability (top 20)")
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(OUT_IMP, dpi=300, bbox_inches="tight")

    report = "\n".join(rep)
    with open(OUT_REPORT, "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nwrote {OUT_RANK}\nwrote {OUT_REPORT}\nwrote {OUT_IMP}")


if __name__ == "__main__":
    main()
