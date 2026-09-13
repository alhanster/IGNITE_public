#!/usr/bin/env python3
"""Builds figure_data/feature_dictionary.csv: for each of the 24 modelled features, records evidence block, condition, definition, coverage, imputation rule, source, observed association with the ranking, and bagging importance.

Runs as part of the `make tables` step, after pu_target_model.py, which writes the importance table joined here. See REPRODUCIBILITY.md, Feature dictionary provenance, for why importances are not recomputed in this script.

Block membership for each feature is imported from attribution_decomposition.py, which partitions the feature matrix at import time.

Inputs:
  data/perturbseq/pu/pu_model_matrix.parquet             coverage (real NaNs) and dtypes
  outputs/model/pu_feature_importance_stability.csv      imp_mean / imp_sd across bags
  figure_data/ranked_atlas.csv                           feature values beside pu_score

Output:
  figure_data/feature_dictionary.csv                     24 rows

Run: PYTHONPATH=src .venv/bin/python src/analysis/tables/build_feature_dictionary.py
"""
import os
import sys

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
from attribution_decomposition import PERTURBATIONAL, GENETIC, OBSERV  # noqa: E402

MATRIX = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
IMPORTANCE = os.path.join(ROOT, "outputs", "model", "pu_feature_importance_stability.csv")
ATLAS = os.path.join(ROOT, "figure_data", "ranked_atlas.csv")
OUT = os.path.join(ROOT, "figure_data", "feature_dictionary.csv")

BLOCK_OF = {}
for _name, _feats in [("genetic", GENETIC), ("observational", OBSERV), ("perturbational", PERTURBATIONAL)]:
    for _f in _feats:
        BLOCK_OF[_f] = _name

# Condition derives from suffix for 12 features; 2 span conditions; zscore_Th* mark Th states, not conditions.
CONDITION_OVERRIDE = {
    "expected_n_regulators_residuals": "across-condition summary (mean of 3)",
    "polar_rank_range": "across-condition summary (range of 3)",
    "has_cytokine": "across-condition summary (coverage)",
}

SOURCE_OF_BLOCK = {
    "genetic": None,  # per-feature; set below
    "observational": "Zhu, Dann et al. 2025 (bioRxiv), Th-vs-Th0 differential expression",
    "perturbational": "Zhu, Dann et al. 2025 (bioRxiv), genome-scale CD4+ T-cell perturb-seq screen",
}
SOURCE_OF_FEATURE = {
    "mis.z_score": "gnomAD v4.1.0 gene constraint metrics",
    "lof.oe_ci.upper": "gnomAD v4.1.0 gene constraint metrics",
    "IEI": "IUIS inborn-errors-of-immunity gene list",
    "gwas_score": "Open Targets Platform 26.06, immune-indication GWAS score",
    "gene_burden_score": "Open Targets Platform 26.06, gene_burden datasource",
}

# See REPRODUCIBILITY.md, Feature dictionary provenance.
DEFINITION = {
    "IEI": "Binary membership of the inborn-errors-of-immunity gene list.",
    "gwas_score": "Per-gene GWAS association score for immune indications. Distinct from the "
                  "broader OT_genetic_immune comparator, and 54-64% zero-inflated by "
                  "drug-status group.",
    "gene_burden_score": "Per-gene rare-variant burden association score for immune "
                         "indications. Distinct from gwas_score, which aggregates "
                         "common-variant credible sets from the same platform release. Fully "
                         "observed but 99.6% zero-inflated: 85 of 19,502 genes carry a "
                         "non-zero score, so a zero records absence of burden evidence "
                         "rather than absence of measurement.",
    "mis.z_score": "Missense constraint Z-score: depletion of missense variation relative to "
                   "expectation.",
    "lof.oe_ci.upper": "LOEUF, the upper bound of the 90% CI on the observed/expected ratio of "
                       "loss-of-function variants. Lower means more constrained.",
    "zscore_Th1": "Differential-expression z-score, Th1-polarized versus resting Th0.",
    "zscore_Th2": "Differential-expression z-score, Th2-polarized versus resting Th0.",
    "zscore_Th17": "Differential-expression z-score, Th17-polarized versus resting Th0.",
    "zscore_Treg": "Differential-expression z-score, Treg-polarized versus resting Th0.",
    "expected_n_regulators_residuals":
        "Trans-regulatory burden. The incoming-regulator count (perturbations significantly "
        "moving the gene's expression) is regressed on log10 mean expression with a Poisson "
        "GLM per condition; the residual is burden in excess of what expression predicts. The "
        "three condition-specific residuals are averaged into one value before modelling.",
    "polar_rank_range": "Range of the three polarization-coefficient ranks across conditions.",
    "has_cytokine": "Coverage indicator: whether any cytokine was measured for the gene. "
                    "Cytokine and cytokine-receptor coverage are coincident in this matrix "
                    "(identical 11,485 genes).",
}
for _c, _lab in [("Rest", "at rest"), ("Stim8hr", "after 8 h stimulation"),
                 ("Stim48hr", "after 48 h stimulation")]:
    DEFINITION[f"polar_coef_rank_{_c}"] = (
        f"Ranked polarization coefficient from the perturbation screen, {_lab}.")
    DEFINITION[f"reg_burden_sig_{_c}"] = (
        f"Flag for a significant lymphocyte-count loss-of-function burden correlation, {_lab}.")
    DEFINITION[f"n_sig_regulated_cytokines_{_c}"] = (
        f"Count of cytokines significantly regulated by perturbation of the gene, {_lab}.")
    DEFINITION[f"n_sig_regulated_cytokine_receptors_{_c}"] = (
        f"Count of cytokine receptors significantly regulated by perturbation of the gene, "
        f"{_lab}.")


def condition_of(feature):
    if feature in CONDITION_OVERRIDE:
        return CONDITION_OVERRIDE[feature]
    for suffix in ("Rest", "Stim8hr", "Stim48hr"):
        if feature.endswith("_" + suffix):
            return suffix
    return ""


def type_of(series):
    """Type comes from values, not dtype: reg_burden_sig_* are float64 flags; IEI/has_cytokine coverage reflects dtype."""
    vals = series.dropna()
    return "binary" if set(vals.unique()) <= {0.0, 1.0} else "continuous"


def main():
    matrix = pd.read_parquet(MATRIX)
    features = GENETIC + OBSERV + PERTURBATIONAL
    missing = [f for f in features if f not in matrix.columns]
    assert not missing, f"features absent from the matrix: {missing}"
    assert len(features) == 24, f"expected 24 features, got {len(features)}"

    n_genes = len(matrix)
    imp = pd.read_csv(IMPORTANCE)
    assert set(imp["feature"]) == set(features), \
        "importance table and the block lists disagree on the feature set"

    atlas = pd.read_csv(ATLAS)

    rows = []
    for f in features:
        col = matrix[f]
        n_ok = int(col.notna().sum())
        block = BLOCK_OF[f]
        # Sign is observed, not model-supplied (gain unsigned); computed only over genes with observed ranking values.
        rho = atlas[[f, "pu_score"]].dropna().corr(method="spearman").iloc[0, 1]
        rows.append({
            "feature": f,
            "block": block,
            "condition": condition_of(f),
            "type": type_of(col),
            "definition": DEFINITION[f],
            "observed_rho": round(float(rho), 4),
            "n_nonmissing": n_ok,
            "pct_nonmissing": round(100.0 * n_ok / n_genes, 2),
            "imputation": ("not imputed (no missing values)" if n_ok == n_genes else
                           "training-row median within each bagging iteration; "
                           "zero if entirely missing within a fit"),
            "source": SOURCE_OF_FEATURE.get(f) or SOURCE_OF_BLOCK[block],
        })

    out = pd.DataFrame(rows).merge(imp, on="feature", how="left", validate="one_to_one")
    assert out["imp_mean"].notna().all(), "importance join left holes"

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}")
    print(f"  {len(out)} features | blocks: {out['block'].value_counts().to_dict()}")
    print(f"  imp_mean sums to {out['imp_mean'].sum():.4f} (normalized gain)")


if __name__ == "__main__":
    main()
