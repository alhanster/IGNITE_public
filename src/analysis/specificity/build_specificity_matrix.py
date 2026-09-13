#!/usr/bin/env python3
"""Computes per-cell bootstrap ROC-AUC with 95 percent confidence intervals for the immune-vs-cardiac specificity control (Fig 7).

For each recovery task (immune-exclusive, cardiac-exclusive held-out trial targets) and each of three scores (PU model, OT genetic_association immune, OT genetic_association cardiac), computes the AUC discriminating that task's exclusive positives from never-labeled genes, using a 2000-draw bootstrap for the confidence interval.

Inputs:
  outputs/model/full_model_pu_scores.csv
  data/specificity/ot_minus_clinical_per_gene.csv
  data/specificity/cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv

Outputs:
  figure_data/specificity_matrix.csv
  figure_data/specificity_panel_n.csv
"""
import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
assert os.path.isfile(os.path.join(REPO_ROOT, "Makefile")), f"repo root not found: {REPO_ROOT}"

sc = pd.read_csv(f"{REPO_ROOT}/outputs/model/full_model_pu_scores.csv")[['gene', 'pu_score', 'pu_role']]
pg = pd.read_csv(f"{REPO_ROOT}/data/specificity/ot_minus_clinical_per_gene.csv")
df = pg.merge(sc, on='gene', how='left').rename(columns={'pu_role': 'immune_role'})

def role(s):
    if s == "Approved":
        return "P"
    if s in ("Phase 1", "Phase 2", "Phase 3"):
        return "trial_heldout"
    return "unlabeled"

cg = pd.read_csv(f"{REPO_ROOT}/data/specificity/cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv")
gcol = 'gene_target' if 'gene_target' in cg.columns else cg.columns[0]
df['cardiac_role'] = df.gene.map(dict(zip(cg[gcol], cg.furthest_stage.map(role)))).fillna("unlabeled")

imm_tt = set(df.loc[df.immune_role == 'trial_heldout', 'gene'])
imm_p = set(df.loc[df.immune_role == 'P', 'gene'])
car_tt = set(df.loc[df.cardiac_role == 'trial_heldout', 'gene'])
car_p = set(df.loc[df.cardiac_role == 'P', 'gene'])
imm_all = imm_tt | imm_p
car_all = car_tt | car_p
ie = imm_tt - car_all
ce = car_tt - imm_all
neither = set(df[(df.immune_role == 'unlabeled') & (df.cardiac_role == 'unlabeled')].gene)
negm = df.gene.isin(neither).values

rng = np.random.default_rng(42)

def cell(scol, posset):
    pm = df.gene.isin(posset).values
    mask = (pm | negm)
    y = pm[mask].astype(int)
    s = df[scol].values[mask]
    keep = ~np.isnan(s)
    y = y[keep]
    s = s[keep]
    a = roc_auc_score(y, s)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    bs = []
    for _ in range(2000):
        idx = np.concatenate([rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)])
        bs.append(roc_auc_score(y[idx], s[idx]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return a, lo, hi, int(y.sum())

rows = []
specs = [("our_model", "pu_score"),
         ("OT_genetic_immune", "immune_genetic_association"),
         ("OT_genetic_cardiac", "cardiac_genetic_association")]

for rec, posset in [("immune_exclusive", ie), ("cardiac_exclusive", ce)]:
    for sname, scol in specs:
        a, lo, hi, n = cell(scol, posset)
        rows.append((rec, sname, round(a, 3), round(lo, 3), round(hi, 3), n,
                     f"{a:.3f} ({lo:.3f}-{hi:.3f})"))

mat = pd.DataFrame(rows, columns=["recover", "score", "auc", "ci_lo", "ci_hi", "n_pos", "ci"])

# rng is one generator consumed in loop order; the removed L2G cells drew last, so the six retained cells' draws are unaffected.

OUT = os.path.join(REPO_ROOT, "figure_data", "specificity_matrix.csv")
mat.to_csv(OUT, index=False)

# max() over per-score n_pos guards against future NaN divergence between scores; currently a no-op.
panel_n = mat.groupby("recover", as_index=False)["n_pos"].max()
PANEL_N = os.path.join(REPO_ROOT, "figure_data", "specificity_panel_n.csv")
panel_n.to_csv(PANEL_N, index=False)

print(f"wrote {os.path.relpath(OUT, REPO_ROOT)}")
print(f"wrote {os.path.relpath(PANEL_N, REPO_ROOT)}")
print(mat.to_string(index=False))
