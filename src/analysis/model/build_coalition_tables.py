#!/usr/bin/env python3
"""Builds the evidence-coalition performance and attribution tables.

Runs after attribution_decomposition.py, reshape_attribution_to_panels.py and compute_nb_significance.py, reading their outputs, reordering columns and dropping one renderer-only column. Writes:

  figure_data/coalition_auc.csv           8 rows, all 2^3 evidence-block coalitions
  figure_data/coalition_ladder.csv        4 rows, model ladder
  figure_data/coalition_significance.csv  3 rows, paired significance tests
  figure_data/coalition_attribution.csv   3 rows, grouped-Shapley shares
  figure_data/coalition_strata.csv        6 rows, per-perturbational-block weak/strong lift
"""
import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"

FIGDATA = os.path.join(ROOT, "figure_data")
COALITIONS = os.path.join(ROOT, "outputs", "model", "functional_genomics_contribution",
                          "attribution_coalition_auc.csv")

SHEETS = {
    "coalition_auc.csv": (COALITIONS, ["coalition", "auc_mean", "auc_sd"]),
    "coalition_ladder.csv": (os.path.join(FIGDATA, "panelA_auc_ladder.csv"),
                             ["model", "label", "auc", "sd", "sem"]),
    "coalition_significance.csv": (os.path.join(FIGDATA, "panelA_significance.csv"),
                                   ["model_a", "model_b", "delta", "nb_p", "stars"]),
    "coalition_attribution.csv": (os.path.join(FIGDATA, "panelB_shapley_shares.csv"),
                                  ["evidence", "shapley_auc", "share_pct"]),
    "coalition_strata.csv": (os.path.join(FIGDATA, "panelE_block_signal.csv"),
                             ["block", "weak_lift", "strong_lift"]),
}


def main():
    out = {}
    for name, (src, cols) in SHEETS.items():
        df = pd.read_csv(src)
        missing = [c for c in cols if c not in df.columns]
        assert not missing, f"{os.path.relpath(src, ROOT)} is missing {missing}"
        # Selecting by name also drops the kind column, used only to tell the renderer which bars to shade.
        out[name] = df[cols]

    coal = out["coalition_auc.csv"]
    assert len(coal) == 8, f"expected all 2^3 coalitions, got {len(coal)}"
    chance = coal.loc[coal["coalition"] == "chance"]
    assert len(chance) == 1 and float(chance["auc_mean"].iloc[0]) == 0.5, \
        "the empty coalition must be present and fixed at chance"

    # Block Shapley values sum exactly to the full model's AUC above chance; mismatch means differing runs.
    attribution = out["coalition_attribution.csv"]
    full_above_chance = float(coal.loc[coal["coalition"] ==
                                       "genetic+observational+perturbational", "auc_mean"].iloc[0]) - 0.5
    total = float(attribution["shapley_auc"].sum())
    assert abs(total - full_above_chance) < 1e-9, \
        f"Shapley values sum to {total}, not the full model's {full_above_chance} above chance"
    assert abs(float(attribution["share_pct"].sum()) - 100.0) < 0.05, "shares do not sum to 100"

    # panelE partitions the perturbational block into six sub-blocks, distinct from the three evidence blocks.
    assert len(out["coalition_strata.csv"]) == 6, "expected 6 perturbational sub-blocks"

    for name, df in out.items():
        df.to_csv(os.path.join(FIGDATA, name), index=False)
        print(f"wrote figure_data/{name}  ({len(df)} rows x {len(df.columns)} cols)")
    print(f"  Shapley efficiency holds: {total:.8f} = full model above chance")


if __name__ == "__main__":
    main()
