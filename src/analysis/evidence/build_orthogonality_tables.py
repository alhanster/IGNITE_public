#!/usr/bin/env python3
"""build_orthogonality_tables.py: correlates functional-genomic features against the five genetic priors, and correlates the genetic priors with each other.

Reads the table produced by build_ranked_atlas_table.py. Performs a join only, with no fitting, seeding, or randomness.

Outputs:
    figure_data/orthogonality_genetic_vs_fg.csv (95 rows: 5 genetic priors x 19 functional-genomic features)
    figure_data/orthogonality_genetic_pairs.csv (10 rows: pairings among the 5 genetic priors)

See REPRODUCIBILITY.md, Other significance tests, for the pairwise-completeness method.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
from attribution_decomposition import PERTURBATIONAL, GENETIC, OBSERV  # noqa: E402

ATLAS = os.path.join(ROOT, "figure_data", "ranked_atlas.csv")
OUT_FG = os.path.join(ROOT, "figure_data", "orthogonality_genetic_vs_fg.csv")
OUT_PAIRS = os.path.join(ROOT, "figure_data", "orthogonality_genetic_pairs.csv")

# Prior versions used for cross-check, kept in drafts/ and untracked; a clean clone lacks them, so the comparison is skipped.
REF_FG = os.path.join(ROOT, "drafts", "genetics_vs_fg_orthogonality.csv")
REF_PAIRS = os.path.join(ROOT, "drafts", "genetic_prior_correlations.csv")

DISPLAY = {
    "lof.oe_ci.upper": "gnomAD LOEUF",
    "mis.z_score": "gnomAD missense z",
    "IEI": "IEI gene membership",
    "gwas_score": "immune-GWAS score",
    "gene_burden_score": "rare-variant burden score",
}

# LOEUF and missense_z index the same depletion-of-variation axis, so their ~50% shared rank variance is expected.
AXIS = {
    "lof.oe_ci.upper": "evolutionary constraint",
    "mis.z_score": "evolutionary constraint",
    "IEI": "rare-variant (curated)",
    "gene_burden_score": "rare-variant (burden)",
    "gwas_score": "common-variant",
}

# Genetic-pairs rows are ordered cross-axis pairs first, the single within-axis pair last.
# The two rare-variant priors are held as SEPARATE axes: curated IEI membership and the
# gene-burden association share 0.34% of rank variance, no more than a cross-axis pair, so
# grouping them would assert a relationship the data does not carry. Contrast LOEUF against
# missense z at ~50%, which is what makes those two one axis rather than two.
PAIR_ORDER = [
    ("lof.oe_ci.upper", "IEI"),
    ("lof.oe_ci.upper", "gene_burden_score"),
    ("lof.oe_ci.upper", "gwas_score"),
    ("mis.z_score", "IEI"),
    ("mis.z_score", "gene_burden_score"),
    ("mis.z_score", "gwas_score"),
    ("IEI", "gene_burden_score"),
    ("IEI", "gwas_score"),
    ("gene_burden_score", "gwas_score"),
    ("lof.oe_ci.upper", "mis.z_score"),
]


def paired(df, a, b):
    """paired holds only rows where both features are observed; every statistic below uses this frame."""
    return df[[a, b]].dropna()


def main():
    atlas = pd.read_csv(ATLAS)
    fg_features = OBSERV + PERTURBATIONAL
    assert len(GENETIC) == 5 and len(fg_features) == 19, \
        f"expected 5 genetic and 19 functional-genomic features, got {len(GENETIC)}/{len(fg_features)}"

    rows = []
    for g in GENETIC:
        for f in fg_features:
            d = paired(atlas, g, f)
            rho, p = spearmanr(d[g], d[f])
            rows.append({"genetic": g, "fg": f, "n": len(d), "rho": rho, "p": p,
                         "shared_var_pct": rho ** 2 * 100})
    fg_tab = pd.DataFrame(rows)
    assert len(fg_tab) == 95, f"expected 95 pairs, got {len(fg_tab)}"

    rows = []
    for a, b in PAIR_ORDER:
        d = paired(atlas, a, b)
        rho, sp = spearmanr(d[a], d[b])
        r, _ = pearsonr(d[a], d[b])
        rows.append({
            "feature_1": DISPLAY[a], "feature_2": DISPLAY[b],
            "axis_1": AXIS[a], "axis_2": AXIS[b],
            "relation": "within-axis" if AXIS[a] == AXIS[b] else "cross-axis",
            "n": len(d),
            # Values are rounded because this sheet is meant to be read as a table, not joined downstream.
            "spearman_rho": round(rho, 4), "spearman_p": f"{sp:.2e}",
            "pearson_r": round(r, 4),
            "shared_variance_pct_spearman": round(rho ** 2 * 100, 2),
        })
    pair_tab = pd.DataFrame(rows)

    _crosscheck(fg_tab, pair_tab)

    fg_tab.to_csv(OUT_FG, index=False)
    pair_tab.to_csv(OUT_PAIRS, index=False)
    print(f"wrote {OUT_FG}\n  {len(fg_tab)} pairs | max |rho| = {fg_tab['rho'].abs().max():.4f} "
          f"| max shared variance = {fg_tab['shared_var_pct'].max():.2f}% "
          f"| pairs above 5%: {(fg_tab['shared_var_pct'] > 5).sum()}")
    print(f"wrote {OUT_PAIRS}\n  {len(pair_tab)} prior pairs | "
          f"{(pair_tab['relation'] == 'within-axis').sum()} within-axis")


def _crosscheck(fg_tab, pair_tab):
    """Cross-checks against drafts/ versions when present; a missing directory skips rather than fails the comparison, and present values must match exactly."""
    if os.path.exists(REF_FG):
        ref = pd.read_csv(REF_FG).sort_values(["genetic", "fg"]).reset_index(drop=True)
        got = fg_tab.sort_values(["genetic", "fg"]).reset_index(drop=True)
        assert list(ref["genetic"]) == list(got["genetic"]), "genetic column disagrees with drafts/"
        assert list(ref["fg"]) == list(got["fg"]), "fg column disagrees with drafts/"
        assert (ref["n"].values == got["n"].values).all(), "n disagrees with drafts/"
        assert np.allclose(ref["rho"], got["rho"], atol=1e-12), "rho disagrees with drafts/"
        print(f"  cross-check: reproduces {os.path.relpath(REF_FG, ROOT)} exactly")

    if os.path.exists(REF_PAIRS):
        ref = pd.read_csv(REF_PAIRS)
        assert len(ref) == len(pair_tab), "prior-pair row count disagrees with drafts/"
        assert (ref["n"].values == pair_tab["n"].values).all(), "prior-pair n disagrees with drafts/"
        assert np.allclose(ref["spearman_rho"], pair_tab["spearman_rho"], atol=5e-5), \
            "prior-pair spearman_rho disagrees with drafts/"
        print(f"  cross-check: reproduces {os.path.relpath(REF_PAIRS, ROOT)} exactly")


if __name__ == "__main__":
    main()
