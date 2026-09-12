"""Backing data for the held-out in-trial ROC curves, full model vs genetics-only, and
the paired DeLong contrast between them.

Computes ROC curves and class counts for the held-out in-trial and non-target genes
and writes them for the renderer as Supplementary Table 7; the class counts are
computed at build time, not hardcoded. Config is imported from the model-analysis
recovery panel, so the genetics-only comparator is re-scored at pum.SEED and
pum.T_BAG. See REPRODUCIBILITY.md, Build-time assertions inside analysis scripts, for the
verification this build performs against panelD_recovery.csv.
"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pu_target_model as pum
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import norm

ROOT    = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUTD    = os.path.join(ROOT, "figure_data")
GENETIC = ["lof.oe_ci.upper", "mis.z_score", "IEI", "gwas_score", "gene_burden_score"]

def _components(y, s):
    pos, neg = s[y == 1], s[y == 0]
    m, n = len(pos), len(neg)
    v01 = np.array([((pos[i] > neg).sum() + 0.5*(pos[i] == neg).sum())/n for i in range(m)])
    v10 = np.array([((pos > neg[j]).sum() + 0.5*(pos == neg[j]).sum())/m for j in range(n)])
    return v01, v10, m, n

def paired_delong(y, s1, s2):
    a1, a2 = roc_auc_score(y, s1), roc_auc_score(y, s2)
    v01_1, v10_1, m, n = _components(y, s1)
    v01_2, v10_2, _, _ = _components(y, s2)
    S = np.cov(np.vstack([v01_1, v01_2]))/m + np.cov(np.vstack([v10_1, v10_2]))/n
    var = S[0,0] + S[1,1] - 2*S[0,1]
    z = (a1 - a2)/np.sqrt(var)
    # 95% CI on the delta
    half = norm.ppf(0.975)*np.sqrt(var)
    return a1, a2, z, 2*(1 - norm.cdf(abs(z))), (a1-a2-half, a1-a2+half)

def auc_ci(y, s):
    a = roc_auc_score(y, s)
    v01, v10, m, n = _components(y, s)
    se = np.sqrt(v01.var(ddof=1)/m + v10.var(ddof=1)/n)
    return a, a - norm.ppf(0.975)*se, a + norm.ppf(0.975)*se

m = pd.read_parquet(pum.MATRIX)
g = pd.read_csv(os.path.join(pum.EVAL, "full_model_pu_scores.csv"))
assert set(g.gene) == set(m.gene) and g.gene.is_unique and m.gene.is_unique
full_sc = m.gene.map(g.set_index("gene").pu_score).values      # committed scores (Fig 4e's file)
assert not np.isnan(full_sc).any()
role  = m["pu_role"].values
P_idx = np.where(role == "P")[0]; U_idx = np.where(role == "unlabeled")[0]
gen_sc, _, _ = pum.pu_bag(m[GENETIC].apply(pd.to_numeric, errors="coerce").values.astype(float),
                          P_idx, U_idx, GENETIC, T=pum.T_BAG, seed=pum.SEED)

sub = (role == "trial_heldout") | (role == "unlabeled")
y   = (role[sub] == "trial_heldout").astype(int)
a_f, a_g, z, p, dci = paired_delong(y, full_sc[sub], gen_sc[sub])
_, f_lo, f_hi = auc_ci(y, full_sc[sub]); _, g_lo, g_hi = auc_ci(y, gen_sc[sub])

os.makedirs(OUTD, exist_ok=True)
rows = []
for name, s in [("full", full_sc[sub]), ("genetics", gen_sc[sub])]:
    fpr, tpr, _ = roc_curve(y, s)
    k = max(1, len(fpr)//2000)                      # thin for a compact CSV
    idx = np.unique(np.r_[np.arange(0, len(fpr), k), len(fpr)-1])
    rows.append(pd.DataFrame({"model": name, "fpr": fpr[idx], "tpr": tpr[idx]}))
pd.concat(rows).to_csv(os.path.join(OUTD, "heldout_trial_roc_curves.csv"), index=False)

# See REPRODUCIBILITY.md, Build-time assertions inside analysis scripts, for the verification performed here.
import make_panelD_recovery as mpd  # noqa: E402

_ev = pd.read_csv(mpd.EVIDENCE)
_imm = mpd.immune_num_by_gene(_ev)
_na = m[m.pu_role != "P"].copy()
_na["sc"] = gen_sc[role != "P"]
_na = _na.sort_values("sc", ascending=False).reset_index(drop=True)
_got = mpd.recovery_counts(_na.gene, _imm)

_panel = pd.read_csv(os.path.join(OUTD, "panelD_recovery.csv"))
_want = dict(zip(_panel["k"].astype(int), _panel["genetics"].astype(int)))

if _got != _want:
    raise SystemExit(
        "genetics-only recovery counts disagree with panelD_recovery.csv\n"
        f"  re-derived here : {_got}\n"
        f"  panelD_recovery : {_want}\n"
        "Both come from a 200-iteration PU-bagging fit on the same five genetic columns at\n"
        "the same seed, so they must agree. A difference means the two re-derivations have\n"
        "diverged -- check pum.SEED / pum.T_BAG and the GENETIC column list in both files."
    )
_verified = "  ".join(f"top-{k}={v}" for k, v in sorted(_got.items()))
print(f"verification OK: genetics-only recovery matches panelD_recovery.csv ({_verified})")

# n_feat_full: non-META column count of the score matrix, as selected in pu_target_model.py.
_n_feat_full = len([c for c in m.columns if c not in pum.META])

meta = {"n_in_trial": int(y.sum()), "n_non_target": int((y == 0).sum()),
        "n_feat_full": _n_feat_full, "n_feat_genetics": len(GENETIC),
        "auc_full": a_f, "auc_full_ci": [f_lo, f_hi],
        "auc_genetics": a_g, "auc_genetics_ci": [g_lo, g_hi],
        "delta": a_f - a_g, "delta_ci95": list(dci),
        "delong_z": z, "delong_p": p,
        "config": "full = committed full_model_pu_scores.csv; genetics-only re-scored at pum.SEED/pum.T_BAG",
        "verification": {
            "claim": "genetics-only ranking reproduces panelD_recovery.csv recovery counts",
            "checked_at_build_time": True,
            "recovery_counts": {str(k): int(v) for k, v in sorted(_got.items())},
        }}
with open(os.path.join(OUTD, "heldout_trial_auc.json"), "w") as fh:
    json.dump(meta, fh, indent=1)
print(json.dumps(meta, indent=1))
