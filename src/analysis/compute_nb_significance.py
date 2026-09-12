#!/usr/bin/env python3
"""See REPRODUCIBILITY.md, Nadeau-Bengio corrected significance (Figure 4 panels a and c)."""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUTD = os.path.join(ROOT, "figure_data")

# Reuses the attribution module's scoring machinery and feature groups to keep results consistent.
sys.path.insert(0, os.path.join(ROOT, "src", "analysis"))
import attribution_decomposition as ad  # noqa: E402

R_SEEDS = ad.SEEDS            # seeds 0..4 (the 5-seed resampling)
K = 5                          # folds per seed (StratifiedKFold(5) in ad.cv_auc)


def stars(p):
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2
            else "*" if p < 5e-2 else "ns")


def nb_pvalue(fold_diffs):
    """Computes the Nadeau-Bengio corrected two-sided p-value and corrected standard error for a paired difference vector."""
    d = np.asarray(fold_diffs, dtype=float)
    n = d.size
    rho = 1.0 / K
    corr = (1.0 / n) + (rho / (1.0 - rho))     # NB variance-inflation factor
    se_nb = np.sqrt(corr * d.var(ddof=1))
    t = d.mean() / se_nb
    p = 2 * stats.t.sf(abs(t), df=n - 1)
    return float(p), float(se_nb)


def main():
    m = pd.read_parquet(ad.MATRIX)
    role = m.pu_role.values
    is_U = role == "unlabeled"
    Pidx = np.where(role == "P")[0]
    Uidx = np.where(is_U)[0]

    # Per-model fold-AUC vectors across the r x k resampling; n = len(SEEDS)*5.
    models = {
        "genetic": ad.GENETIC,
        "gen+obs": ad.GENETIC + ad.OBSERV,
        "gen+perturbational": ad.GENETIC + ad.PERTURBATIONAL,
        "full": ad.GENETIC + ad.PERTURBATIONAL + ad.OBSERV,
    }
    fa = {name: [] for name in models}
    for s in R_SEEDS:
        for name, cols in models.items():
            fa[name].extend(ad.cv_auc(m, cols, Pidx, Uidx, is_U, seed=s))
    fa = {k: np.array(v) for k, v in fa.items()}

    diffs = {
        "gen+obs": fa["gen+obs"] - fa["genetic"],
        "gen+perturbational": fa["gen+perturbational"] - fa["genetic"],
        "full": fa["full"] - fa["genetic"],
    }
    nb = {}
    for k, d in diffs.items():
        p, se = nb_pvalue(d)
        nb[k] = {"delta": float(d.mean()), "p": p, "se": se, "stars": stars(p)}

    # Patches panelA_significance.csv: keeps reshape's delta, adds NB p-value and significance stars.
    pa_path = os.path.join(OUTD, "panelA_significance.csv")
    pa = pd.read_csv(pa_path)
    key_by_b = {"genetic+observational": "gen+obs",
                "genetic+perturbational": "gen+perturbational", "full": "full"}
    # Writes the p-value as nb_p rather than overwriting in place, so the file cannot be mistaken for Wilcoxon p-values.
    pa["nb_p"] = pa["model_b"].map(lambda b: nb[key_by_b[b]]["p"])
    pa["stars"] = pa["model_b"].map(lambda b: nb[key_by_b[b]]["stars"])
    pa = pa.drop(columns=[c for c in ("wilcoxon_p",) if c in pa.columns])
    pa.to_csv(pa_path, index=False)

    # Patches panelC_marginal.csv: keeps reshape's delta, adds NB standard error and p-value.
    pc_path = os.path.join(OUTD, "panelC_marginal.csv")
    pc = pd.read_csv(pc_path)
    ev = {"perturbational": "gen+perturbational", "observational": "gen+obs", "full": "full"}
    pc["sem"] = pc["evidence"].map(lambda e: nb[ev[e]]["se"])
    pc["nb_p"] = pc["evidence"].map(lambda e: nb[ev[e]]["p"])
    pc = pc.drop(columns=[c for c in ("wilcoxon_p",) if c in pc.columns])
    pc.to_csv(pc_path, index=False)

    print("Nadeau-Bengio corrected significance written (n=%d folds, rho=1/%d):"
          % (len(next(iter(fa.values()))), K))
    for name, r in nb.items():
        print("  %-11s delta=%+.4f  p=%.4f  %s" % (name, r["delta"], r["p"], r["stars"]))
    print("patched: panelA_significance.csv, panelC_marginal.csv")


if __name__ == "__main__":
    main()
