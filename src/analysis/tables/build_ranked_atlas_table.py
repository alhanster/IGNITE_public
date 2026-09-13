#!/usr/bin/env python3
"""Assembles the per-gene ranked atlas: PU score, rank percentile, drug-status label, and
the 24 PU-model feature values, joined on gene, from outputs/model/full_model_pu_scores.csv
and data/perturbseq/pu/pu_model_matrix.parquet. Writes figure_data/ranked_atlas.csv.

Run: PYTHONPATH=src .venv/bin/python src/analysis/tables/build_ranked_atlas_table.py

See REPRODUCIBILITY.md, Stage-1 step order, for build-order requirements.
"""
import os
import sys
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
import pu_target_model as pum  # owns META; not reimplemented (see build_heldout_trial_auc.py)

SCORES = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
MATRIX = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
OUT = os.path.join(ROOT, "figure_data", "ranked_atlas.csv")

# pu_role/furthest_stage deduped to score table copy; crossdonor_correlation_mean kept as own meta column.
DUPLICATE_META = ["pu_role", "furthest_stage"]


def drug_status(row):
    """Canonical drug-status partition, matches build_score_by_group.py::assign_group."""
    if row["pu_role"] == "P":
        return "approved"
    return "in-trial" if str(row["furthest_stage"]).startswith("Phase") else "non-target"


def main():
    scores = pd.read_csv(SCORES)
    matrix = pd.read_parquet(MATRIX).drop(columns=DUPLICATE_META)

    missing = {"gene", "pu_score", "rank_pctile", "rank", "pu_role", "furthest_stage",
               "crossdonor_confidence"} - set(scores.columns)
    if missing:
        raise SystemExit(f"full_model_pu_scores.csv missing columns: {sorted(missing)}")
    assert scores["gene"].is_unique and matrix["gene"].is_unique
    assert set(scores["gene"]) == set(matrix["gene"]), "score table and feature matrix disagree on gene set"

    out = scores.merge(matrix, on="gene", how="left", validate="one_to_one")
    out["drug_status"] = out.apply(drug_status, axis=1)

    feature_cols = [c for c in matrix.columns if c not in pum.META]
    assert len(feature_cols) == 24, f"expected 24 model features, got {len(feature_cols)}"
    cols = (["gene", "pu_score", "rank_pctile", "rank", "drug_status", "pu_role",
             "furthest_stage", "crossdonor_confidence", "crossdonor_correlation_mean"]
            + feature_cols)
    out = out[cols].sort_values("rank").reset_index(drop=True)

    # Published under a name that says which indications the stage is for: data/specificity/ carries a
    # cardiac furthest_stage, and the Open Targets inputs this reads keep the plain name internally.
    out = out.rename(columns={"furthest_stage": "furthest_imm_stage"})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}")
    print(f"  {len(out)} genes, {len(feature_cols)} model features, "
          f"drug_status counts: {out['drug_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
