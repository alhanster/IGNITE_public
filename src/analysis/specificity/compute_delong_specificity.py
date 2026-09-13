#!/usr/bin/env python3
"""DeLong paired-ROC significance tests for the immune-vs-cardiac specificity control (Fig 7), computed from on-disk per-gene data with no network calls.

DeLong's test compares two ROC curves over the same cases, so it applies only within a recovery panel, where the two scores share identical positives and negatives. Cross-panel comparisons use a bootstrap of the difference instead.

delong_test's signature is also used by build_leakage_controlled.py.

Inputs (repo-relative):
  outputs/model/full_model_pu_scores.csv                                     immune roles and pu_score
  data/specificity/ot_minus_clinical_per_gene.csv                            OT genetic-association scores
  data/specificity/cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv    cardiac labels
Output:
  figure_data/delong_specificity_tests.csv (panel, A, B, auc_A, auc_B, delta, z, p, p_holm, sig)

Method references: DeLong et al. 1988; fast midrank covariance from Sun and Xu 2014.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from _stats import holm  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))

PU_PATH      = os.path.join(REPO_ROOT, "outputs", "model", "full_model_pu_scores.csv")
OT_PATH      = os.path.join(REPO_ROOT, "data", "specificity",
                            "ot_minus_clinical_per_gene.csv")
CARDIAC_PATH = os.path.join(REPO_ROOT, "data", "specificity",
                            "cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv")
OUT_PATH     = os.path.join(REPO_ROOT, "figure_data", "delong_specificity_tests.csv")


# ---------------------------------------------------------------------------
# Fast DeLong (Sun and Xu 2014)
# ---------------------------------------------------------------------------
def _compute_midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x)
    T = np.zeros(N, dtype=float); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float); T2[J] = T
    return T2


def _fastdelong(preds_pos, preds_neg):
    """preds_* : (2, n) arrays; rows are the two predictors' scores."""
    m = preds_pos.shape[1]; n = preds_neg.shape[1]; k = preds_pos.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r] = _compute_midrank(preds_pos[r])
        ty[r] = _compute_midrank(preds_neg[r])
        tz[r] = _compute_midrank(np.r_[preds_pos[r], preds_neg[r]])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    return aucs, np.atleast_2d(cov)


def delong_test(score_a, score_b, y):
    pos = y == 1; neg = y == 0
    pp = np.vstack([score_a[pos], score_b[pos]])
    pn = np.vstack([score_a[neg], score_b[neg]])
    aucs, cov = _fastdelong(pp, pn)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    z = (aucs[0] - aucs[1]) / np.sqrt(var) if var > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z))
    return aucs, z, p


def stars(p):
    # Uses the Holm-corrected p; fig_specificity.R's legend must list every symbol drawn, including **.
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 0.05: return "*"
    return "ns"


def role_from_stage(s):
    if s == "Approved": return "P"
    if s in ("Phase 1", "Phase 2", "Phase 3"): return "trial_heldout"
    return "unlabeled"   # Unknown-stage -> non-target, matching the immune build


def main():
    for pth in (PU_PATH, OT_PATH, CARDIAC_PATH):
        if not os.path.exists(pth):
            raise FileNotFoundError(pth)

    pu   = pd.read_csv(PU_PATH)[["gene", "pu_score", "pu_role"]]
    otpg = pd.read_csv(OT_PATH)[["gene", "immune_genetic_association",
                                 "cardiac_genetic_association"]]
    card = pd.read_csv(CARDIAC_PATH)
    card_role = {g: role_from_stage(s)
                 for g, s in zip(card["gene_target"], card["furthest_stage"])}

    df = pu.merge(otpg, on="gene", how="left")
    df["immune_role"]  = df["pu_role"]
    df["cardiac_role"] = df["gene"].map(card_role).fillna("unlabeled")
    for c in ("immune_genetic_association", "cardiac_genetic_association"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    imm_tt = set(df.loc[df.immune_role == "trial_heldout", "gene"])
    imm_p  = set(df.loc[df.immune_role == "P", "gene"])
    car_tt = set(df.loc[df.cardiac_role == "trial_heldout", "gene"])
    car_p  = set(df.loc[df.cardiac_role == "P", "gene"])
    imm_excl = imm_tt - (car_tt | car_p)
    car_excl = car_tt - (imm_tt | imm_p)
    neither  = set(df[(df.immune_role == "unlabeled") &
                      (df.cardiac_role == "unlabeled")].gene)
    neg_mask = df.gene.isin(neither).values

    print(f"immune-exclusive : {len(imm_excl)}")
    print(f"cardiac-exclusive: {len(car_excl)}")
    print(f"clean negatives  : {int(neg_mask.sum())}")

    def vec(posset, scol):
        pm = df.gene.isin(posset).values
        mask = pm | neg_mask
        return df[scol].values[mask].astype(float), pm[mask].astype(int)

    specs = {"our_model": "pu_score",
             "OT_genetic_immune": "immune_genetic_association",
             "OT_genetic_cardiac": "cardiac_genetic_association"}
    pairs = [("our_model", "OT_genetic_immune"),
             ("our_model", "OT_genetic_cardiac"),
             ("OT_genetic_immune", "OT_genetic_cardiac")]

    rows = []
    for panel, posset in [("immune_exclusive", imm_excl),
                          ("cardiac_exclusive", car_excl)]:
        _, y = vec(posset, "pu_score")
        # See REPRODUCIBILITY.md, Multiple-testing correction.
        raw = []
        for a, b in pairs:
            sa, _ = vec(posset, specs[a])
            sb, _ = vec(posset, specs[b])
            raw.append((a, b) + delong_test(sa, sb, y))
        adj = holm([r[4] for r in raw])
        for (a, b, aucs, z, p), p_h in zip(raw, adj):
            # Correction is applied to the raw float p-values here, before rounding to two significant figures for the CSV.
            rows.append({"panel": panel, "A": a, "B": b,
                         "auc_A": round(aucs[0], 3), "auc_B": round(aucs[1], 3),
                         "delta": round(aucs[0] - aucs[1], 3),
                         "z": round(z, 2), "p": f"{p:.2e}", "p_holm": f"{p_h:.2e}",
                         "sig": stars(p_h)})

    out = pd.DataFrame(rows)
    out.to_csv(OUT_PATH, index=False)
    print(out.to_string(index=False))
    print("Saved:", OUT_PATH)


if __name__ == "__main__":
    main()
