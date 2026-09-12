"""Label-permutation null helpers; ladder from figure_data/panelA_auc_ladder.csv, not the coalition csv."""
import os
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FIGDATA = os.path.join(ROOT, "figure_data")
LADDER_CSV = os.path.join(FIGDATA, "panelA_auc_ladder.csv")
NULL_CSV = os.path.join(FIGDATA, "label_permutation_null.csv")
PVAL_CSV = os.path.join(FIGDATA, "label_permutation_pvalues.csv")

# Checkpoints go to gitignored scratch, not figure_data/: verify-tables fails on any untracked file there.
CKPT_DIR = os.path.join(ROOT, "outputs", "permutation")

RUNGS = ["genetic", "genetic+observational", "genetic+perturbational", "full"]

# Maps panelA_auc_ladder.csv's display labels to the internal rung names.
_LADDER_KEY = {"genetic only": "genetic", "+ observational": "genetic+observational",
               "+ perturbational": "genetic+perturbational", "full": "full"}


def observed_ladder():
    """The four observed CV AUCs, keyed by rung name."""
    if not os.path.exists(LADDER_CSV):
        raise FileNotFoundError(
            f"{LADDER_CSV}\nThe observed ladder is committed; restore it from git or run "
            "`make tables` to rewrite it.")
    df = pd.read_csv(LADDER_CSV)
    obs = {_LADDER_KEY[m]: float(a) for m, a in zip(df["model"], df["auc"])
           if m in _LADDER_KEY}
    missing = set(RUNGS) - set(obs)
    if missing:
        raise SystemExit(f"panelA_auc_ladder.csv is missing rungs: {sorted(missing)}")
    return obs


def pvalue_table(null, obs):
    """One-sided empirical p per rung, relative to the observed ladder.

    p = (#{null >= observed} + 1) / (n + 1); with 1,000 permutations the smallest reportable p is 0.000999, never an exact 0.

    Levels only. Permuting the labels drives every rung to chance, so a difference between two rungs
    tested against this null asks whether an increment could arise when nothing predicts anything,
    not whether one rung beats another; scrambled_feature_control.py is that test. Increment rows
    were reported here until they were removed as not a valid comparison -- do not re-add them.
    """
    rows, n = [], len(null)
    def row(metric, kind, o, nd):
        return dict(metric=metric, kind=kind, observed=o,
                    null_mean=float(nd.mean()), null_sd=float(nd.std()),
                    null_p95=float(np.percentile(nd, 95)),
                    emp_p=(np.sum(nd >= o) + 1) / (n + 1), n_perm=n)
    for name in RUNGS:
        rows.append(row(name, "level", obs[name], null[name].values))
    return pd.DataFrame(rows)
