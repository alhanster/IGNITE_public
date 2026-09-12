#!/usr/bin/env python
"""Decomposes PU model performance across three feature groups, genetic, perturbational, and observational, using grouped Shapley values, paired significance tests, and held-out-target recovery.

The feature groups form an exact partition of the model features, defined in the GENETIC, PERTURBATIONAL, and OBSERV lists below. GENETIC covers evolutionary and population priors. PERTURBATIONAL covers perturbation-derived effects (CRISPRi knockdown). OBSERVATIONAL covers a gene's own expression change during polarization.

Computes grouped Shapley values over AUC-above-chance across all 8 coalitions of the three groups (5-fold PU-bagging cross-validation, averaged over 5 seeds), paired Wilcoxon signed-rank tests on matched per-fold AUCs, and top-k recovery of held-out trial targets.

See REPRODUCIBILITY.md, Nadeau-Bengio corrected significance (Figure 4 panels a and c), for the corrected significance values used in published figures.

Input: data/perturbseq/pu/pu_model_matrix.parquet
Outputs: outputs/model/functional_genomics_contribution/attribution_results.json, attribution_coalition_auc.csv, outputs/model/figure_attribution.png
"""
import os
import json
from itertools import combinations
from math import factorial
import numpy as np
import pandas as pd
from scipy import stats
from xgboost import XGBClassifier

# Aborts if xgboost version differs, since mismatches silently alter gene-level results.
import sys as _sys, os as _os
_r = _os.path.abspath(__file__)
while _r != _os.path.dirname(_r) and not _os.path.isfile(_os.path.join(_r, "Makefile")):
    _r = _os.path.dirname(_r)
_sys.path.insert(0, _os.path.join(_r, "src", "analysis"))
from _version_guard import check_pins  # noqa: E402
check_pins()
# -----------------------------------------------------------------------------
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
PU_DIR = os.path.join(ROOT, "data", "perturbseq", "pu")
OUT = os.path.join(ROOT, "outputs", "model", "functional_genomics_contribution")
os.makedirs(OUT, exist_ok=True)
MATRIX = os.path.join(PU_DIR, "pu_model_matrix.parquet")
OUT_FIG = os.path.join(ROOT, "outputs", "model", "figure_attribution.png")

META = ["gene", "pu_role", "furthest_stage", "crossdonor_correlation_mean"]
GENETIC = ["lof.oe_ci.upper", "mis.z_score", "IEI", "gwas_score", "gene_burden_score"]  # gwas_score_height dropped (non-immune height GWAS proxy)
# Lists must exactly partition non-META columns: an omitted feature drops silently from Shapley, not from production.
PERTURBATIONAL = ["expected_n_regulators_residuals",
          "polar_coef_rank_Rest", "polar_coef_rank_Stim8hr", "polar_coef_rank_Stim48hr",
          "polar_rank_range",
          "n_sig_regulated_cytokines_Rest", "n_sig_regulated_cytokines_Stim8hr",
          "n_sig_regulated_cytokines_Stim48hr", "has_cytokine",
          "n_sig_regulated_cytokine_receptors_Rest", "n_sig_regulated_cytokine_receptors_Stim8hr",
          "n_sig_regulated_cytokine_receptors_Stim48hr",
          "reg_burden_sig_Rest", "reg_burden_sig_Stim8hr", "reg_burden_sig_Stim48hr"]
OBSERV = ["zscore_Th1", "zscore_Th2", "zscore_Th17", "zscore_Treg"]
GROUPS = {"genetic": GENETIC, "perturbational": PERTURBATIONAL, "observational": OBSERV}
G = list(GROUPS.keys())
# See REPRODUCIBILITY.md, Cross-validation and seeds.
SEEDS = [0, 1, 2, 3, 4]
N_SEED = len(SEEDS)
T_BAG = 30


_MAX_DEPTH = int(os.environ.get("MODEL_MAX_DEPTH", "2"))   # final config = 2 (2026-07-23)


def base_learner(seed):
    # max_depth=2 matches the production base_learner config; MODEL_MAX_DEPTH env var overrides it for depth-diff runs.
    return XGBClassifier(n_estimators=120, max_depth=_MAX_DEPTH, learning_rate=0.1, subsample=0.8,
                         colsample_bytree=0.8, eval_metric="logloss", n_jobs=4,
                         verbosity=0, random_state=seed)


def cv_auc(m, cols, Pidx, Uidx, is_U, seed):
    if len(cols) == 0:
        return [0.5] * 5
    X = m[cols].apply(pd.to_numeric, errors="coerce").values.astype(float)
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    out = []
    for trp, tep in skf.split(Pidx, np.ones(len(Pidx))):
        trainP, testP = Pidx[trp], Pidx[tep]
        rr = np.random.RandomState(seed * 13 + 7)
        ssum = np.zeros(len(m)); scnt = np.zeros(len(m))
        for t in range(T_BAG):
            neg = rr.choice(Uidx, len(trainP), replace=True)
            tri = np.concatenate([trainP, neg]); yt = np.r_[np.ones(len(trainP)), np.zeros(len(neg))]
            md = np.nanmedian(X[tri], axis=0); md = np.where(np.isnan(md), 0, md)
            Xtr = np.where(np.isnan(X[tri]), md, X[tri]); Xall = np.where(np.isnan(X), md, X)
            c = base_learner(t).fit(Xtr, yt); p = c.predict_proba(Xall)[:, 1]
            oob = np.ones(len(m), bool); oob[neg] = False; ssum[oob] += p[oob]; scnt[oob] += 1
        s = np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan)
        msk = np.zeros(len(m), bool); msk[testP] = True; msk[is_U] = True
        y = np.zeros(len(m)); y[testP] = 1; v = msk & ~np.isnan(s)
        out.append(roc_auc_score(y[v], s[v]))
    return out


def full_score(m, cols, Pidx, Uidx, seed=0, T=100):
    if len(cols) == 0:
        return np.full(len(m), np.nan)
    X = m[cols].apply(pd.to_numeric, errors="coerce").values.astype(float)
    rng = np.random.RandomState(seed); ssum = np.zeros(len(m)); scnt = np.zeros(len(m))
    for t in range(T):
        neg = rng.choice(Uidx, len(Pidx), replace=True)
        tri = np.concatenate([Pidx, neg]); yt = np.r_[np.ones(len(Pidx)), np.zeros(len(neg))]
        md = np.nanmedian(X[tri], axis=0); md = np.where(np.isnan(md), 0, md)
        Xtr = np.where(np.isnan(X[tri]), md, X[tri]); Xall = np.where(np.isnan(X), md, X)
        c = base_learner(t).fit(Xtr, yt); p = c.predict_proba(Xall)[:, 1]
        oob = np.ones(len(m), bool); oob[neg] = False; ssum[oob] += p[oob]; scnt[oob] += 1
    return np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan)


def main():
    m = pd.read_parquet(MATRIX)
    allf = [c for c in m.columns if c not in META]
    assert set(GENETIC + PERTURBATIONAL + OBSERV) == set(allf), "groups must partition all features"
    role = m.pu_role.values
    is_P, is_U, is_T = role == "P", role == "unlabeled", role == "trial_heldout"
    Pidx, Uidx = np.where(is_P)[0], np.where(is_U)[0]

    def cols_for(coal):
        # Iterates G, not the frozenset coalitions: hash-randomized order would change colsample_bytree=0.8 sampling and AUCs.
        c = []
        for g in G:
            if g in coal:
                c += GROUPS[g]
        return c

    coalitions = [frozenset(c) for k in range(len(G) + 1) for c in combinations(G, k)]
    coal_auc = {}
    for coal in coalitions:
        aucs = []
        for s in SEEDS:
            aucs += cv_auc(m, cols_for(coal), Pidx, Uidx, is_U, seed=s)
        coal_auc[coal] = np.array(aucs)

    def v(coal):
        return coal_auc[frozenset(coal)].mean() - 0.5

    n = len(G)
    shap = {}
    for g in G:
        others = [x for x in G if x != g]
        tot = 0.0
        for k in range(len(others) + 1):
            for S in combinations(others, k):
                S = set(S)
                w = factorial(len(S)) * factorial(n - len(S) - 1) / factorial(n)
                tot += w * (v(S | {g}) - v(S))
        shap[g] = tot
    V = v(set(G))

    full = coal_auc[frozenset(G)]
    prior = coal_auc[frozenset({"genetic"})]
    prior_obs = coal_auc[frozenset({"genetic", "observational"})]
    prior_perturbational = coal_auc[frozenset({"genetic", "perturbational"})]

    def wp(a, b):
        return float(stats.wilcoxon(a, b)[1])

    # held-out recovery at top-k
    scf = np.nanmean([full_score(m, cols_for(G), Pidx, Uidx, seed=s) for s in SEEDS], axis=0)
    scp = np.nanmean([full_score(m, GENETIC, Pidx, Uidx, seed=s) for s in SEEDS], axis=0)
    scpo = np.nanmean([full_score(m, GENETIC + OBSERV, Pidx, Uidx, seed=s) for s in SEEDS], axis=0)

    def hits(score, ks=(50, 100, 200, 500)):
        pool = is_T | is_U; idx = np.where(pool)[0]
        order = idx[np.argsort(-np.nan_to_num(score[idx], nan=-1))]
        return {k: int(is_T[order[:k]].sum()) for k in ks}
    hf, hp, hpo = hits(scf), hits(scp), hits(scpo)

    results = {
        "groups": {k: len(v_) for k, v_ in GROUPS.items()},
        # n_auc_samples = len(SEEDS) * n_folds; panelA sem = sd / sqrt(n_auc_samples).
        "n_auc_samples": int(len(coal_auc[frozenset(G)])),
        "coalition_auc": {"+".join(sorted(c)) or "chance": {"mean": float(coal_auc[c].mean()),
                          "sd": float(coal_auc[c].std())} for c in coalitions},
        "full_above_chance": float(V),
        "shapley_above_chance": {g: float(shap[g]) for g in G},
        "shapley_share_pct": {g: float(100 * shap[g] / V) for g in G},
        "paired_deltas": {
            # panelC sem = population std (ddof=0) of per-fold paired deltas / sqrt(n_folds), matching panelA.
            "full_vs_genetic": {"delta": float((full - prior).mean()),
                                "sem": float((full - prior).std() / np.sqrt(len(full))),
                                "wilcoxon_p": wp(full, prior)},
            "perturbational_over_genetic": {"delta": float((prior_perturbational - prior).mean()),
                                    "sem": float((prior_perturbational - prior).std() / np.sqrt(len(prior_perturbational))),
                                    "wilcoxon_p": wp(prior_perturbational, prior)},
            "observational_over_genetic": {"delta": float((prior_obs - prior).mean()),
                                           "sem": float((prior_obs - prior).std() / np.sqrt(len(prior_obs))),
                                           "wilcoxon_p": wp(prior_obs, prior)},
        },
        "trial_recovery_topk": [{"k": k, "genetic": hp[k], "genetic_obs": hpo[k], "full": hf[k],
                                 "gain_full_vs_genetic": hf[k] - hp[k]} for k in (50, 100, 200, 500)],
    }
    json.dump(results, open(os.path.join(OUT, "attribution_results.json"), "w"), indent=2)
    pd.DataFrame([{"coalition": "+".join(sorted(c)) or "chance", "auc_mean": coal_auc[c].mean(),
                   "auc_sd": coal_auc[c].std()} for c in coalitions]).to_csv(
        os.path.join(OUT, "attribution_coalition_auc.csv"), index=False)
    make_figure(results, coal_auc, shap, V, full, prior, prior_obs, prior_perturbational)
    print(json.dumps({"shapley_share_pct": results["shapley_share_pct"],
                      "full_above_chance": V,
                      "paired": results["paired_deltas"]}, indent=1))


def make_figure(results, coal_auc, shap, V, full, prior, prior_obs, prior_perturbational):
    try:
        import sys
        sys.path.insert(0, os.path.expanduser("~"))
        from figure_style_kernel import apply_figure_style, goodness_arrow  # if packaged
    except Exception:
        apply_figure_style = None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if apply_figure_style:
        apply_figure_style(sizes=(8, 7, 6))
    C = {"genetic": "#999999", "perturbational": "#009E73", "observational": "#0072B2"}
    order = ["genetic", "perturbational", "observational"]
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2))

    ax = axes[0]; bottom = 0
    for g in order:
        val = shap[g]; ax.bar(0, val, bottom=bottom, color=C[g], width=0.6, edgecolor="white")
        ax.text(0, bottom + val / 2, f"{g}\n{100*shap[g]/V:.0f}%", ha="center", va="center",
                fontsize=8, color="white" if g != "genetic" else "black", fontweight="bold")
        bottom += val
    ax.set_xlim(-0.6, 0.6); ax.set_xticks([]); ax.set_ylabel("Shapley value (AUC above chance)")
    ax.set_title("Share of model performance\nby evidence type")
    ax.text(0, bottom + 0.008, f"total {V:.3f}", ha="center", fontsize=7)

    ax = axes[1]
    seq = [("genetic",), ("genetic", "observational"), ("genetic", "perturbational", "observational")]
    labels = ["genetic\nonly", "+ observational", "+ perturbational\n(full)"]
    means = [coal_auc[frozenset(c)].mean() for c in seq]
    n_obs = len(full)
    sems = [coal_auc[frozenset(c)].std() / np.sqrt(n_obs) for c in seq]
    xs = np.arange(3)
    ax.errorbar(xs, means, yerr=sems, fmt="o-", color="#333333", capsize=3, markersize=6, lw=1.2)
    for x, mn in zip(xs, means):
        ax.text(x, mn + 0.0012, f"{mn:.3f}", ha="center", fontsize=7)
    ax.set_xticks(xs); ax.set_xticklabels(labels); ax.set_ylabel("CV ROC-AUC (approved vs unlabeled)")
    ax.set_ylim(0.760, 0.792); ax.set_title("Adding functional genomics\nto genetic priors")
    ax.text(0.03, 0.04, "error = SEM across folds;\ncomparison is paired (see panel 3)",
            transform=ax.transAxes, fontsize=6, color="0.4", va="bottom")

    ax = axes[2]; xmax = 0.030
    contribs = [("genetic+perturbational\nvs genetic", prior_perturbational - prior, C["perturbational"], prior_perturbational),
                ("genetic+obs\nvs genetic", prior_obs - prior, C["observational"], prior_obs),
                ("full\nvs genetic", full - prior, "#333333", full)]
    yy = np.arange(len(contribs))[::-1]
    for y, (lab, d, col, arm) in zip(yy, contribs):
        ax.barh(y, d.mean(), xerr=d.std() / np.sqrt(n_obs), color=col, capsize=3, height=0.55)
        p = float(stats.wilcoxon(arm, prior)[1])
        ax.text(xmax * 0.98, y, f"p={p:.0e}".replace("e-0", "e-"), va="center", ha="right",
                fontsize=6.5, color="0.3")
    ax.set_yticks(yy); ax.set_yticklabels([c[0] for c in contribs])
    ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("ΔAUC over genetic-only (SEM)")
    ax.set_title("Marginal contribution\n(paired across folds)"); ax.set_xlim(0, xmax)

    fig.suptitle("How much does perturbational functional genomics add to genetic priors?", fontsize=12, y=1.02)
    fig.tight_layout(); fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    print("wrote", OUT_FIG)


if __name__ == "__main__":
    main()
