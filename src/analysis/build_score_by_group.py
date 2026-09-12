#!/usr/bin/env python3
"""build_score_by_group.py

Derives the backing table for Figure 4 panel d directly from the PU model output.

Group definition: pu_role == 'P' -> approved; furthest_stage starting with 'Phase' -> in-trial; else -> non-target. See REPRODUCIBILITY.md, Drug status table, for the redundant-key check.

Input : outputs/model/full_model_pu_scores.csv (produced by src/analysis/pu_target_model.py)
Output: figure_data/score_by_group.csv, figure_data/score_by_group_stats.json

Run: python src/analysis/build_score_by_group.py
"""
import json
import os
import sys

import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _stats import holm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SCORES = os.path.join(REPO, "outputs", "model", "full_model_pu_scores.csv")
OUT_DIR = os.path.join(REPO, "figure_data")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_CSV = os.path.join(OUT_DIR, "score_by_group.csv")
OUT_JSON = os.path.join(OUT_DIR, "score_by_group_stats.json")

GROUP_ORDER = ["non-target", "in-trial", "approved"]


def assign_group(row):
    if row.pu_role == "P":
        return "approved"
    return "in-trial" if str(row.furthest_stage).startswith("Phase") else "non-target"


def main():
    gp = pd.read_csv(SCORES)
    missing = {"gene", "pu_score", "rank_pctile", "pu_score_cv", "rank_pctile_cv",
               "pu_role", "furthest_stage"} - set(gp.columns)
    if missing:
        raise SystemExit(f"full_model_pu_scores.csv missing columns: {sorted(missing)}")

    gp["grp"] = gp.apply(assign_group, axis=1)
    out = gp[["gene", "grp", "rank_pctile", "pu_score", "rank_pctile_cv", "pu_score_cv"]].copy()
    # See REPRODUCIBILITY.md, What a green gate does not cover, for the OOB filtering versus test-frame note.
    keep = out["rank_pctile_cv"].notna()
    # The tests and the JSON's group counts run on the unfiltered frame, so the two scores must
    # cover the same genes or panel d's axis labels overstate what is plotted.
    if not keep.equals(out["rank_pctile"].notna()):
        raise SystemExit("rank_pctile and rank_pctile_cv disagree on which genes are scored; "
                         "the published and cross-fitted scores must cover the same genes")
    out[keep].to_csv(OUT_CSV, index=False)

    groups = [out.loc[out.grp == g, "pu_score_cv"].values for g in GROUP_ORDER]
    H, p_kw = stats.kruskal(*groups)
    F, p_an = stats.f_oneway(*groups)

    # All three groups are scored out of sample on the cross-fitted score; see REPRODUCIBILITY.md, Other significance tests.
    U_ho, p_ho = stats.mannwhitneyu(
        out.loc[out.grp == "in-trial", "pu_score_cv"],
        out.loc[out.grp == "non-target", "pu_score_cv"],
        alternative="greater",
    )

    stats_blob = {
        "heldout_in_trial_vs_non_target": {
            "U": float(U_ho),
            "p_one_sided": float(p_ho),
            "note": "one-sided companion of the in-trial vs non-target pairwise test",
        },
        "score_basis": "cross-fitted, 5 folds over positives at T=200",
        "groups": {
            g: {
                "n": int((out.grp == g).sum()),
                "median": float(out.loc[out.grp == g, "rank_pctile_cv"].median()),
                "mean": float(out.loc[out.grp == g, "rank_pctile_cv"].mean()),
                "q25": float(out.loc[out.grp == g, "rank_pctile_cv"].quantile(0.25)),
                "q75": float(out.loc[out.grp == g, "rank_pctile_cv"].quantile(0.75)),
            }
            for g in GROUP_ORDER
        },
        "kruskal_H": float(H),
        "kruskal_p": float(p_kw),
        "anova_F": float(F),
        "anova_p": float(p_an),
        "pairwise_mannwhitney": [],
        "_provenance": "derived from outputs/model/full_model_pu_scores.csv by build_score_by_group.py",
    }

    pairs = []
    for i in range(len(GROUP_ORDER)):
        for j in range(i + 1, len(GROUP_ORDER)):
            a, b = GROUP_ORDER[i], GROUP_ORDER[j]
            U, p = stats.mannwhitneyu(
                out.loc[out.grp == a, "pu_score_cv"], out.loc[out.grp == b, "pu_score_cv"],
                alternative="two-sided",
            )
            pairs.append((a, b, float(U), float(p)))

    # Holm applies only to pairwise tests (kruskal_p, anova_p uncorrected); column is p_holm not p_bonf.
    for (a, b, U, p), p_holm in zip(pairs, holm([q[3] for q in pairs])):
        # The significance ladder applies to the corrected p-value.
        stars = ("****" if p_holm < 1e-4 else "***" if p_holm < 1e-3
                 else "**" if p_holm < 1e-2 else "*" if p_holm < 0.05 else "ns")
        stats_blob["pairwise_mannwhitney"].append(
            {"a": a, "b": b, "U": U, "p": p, "p_holm": p_holm, "stars": stars}
        )

    with open(OUT_JSON, "w") as fh:
        json.dump(stats_blob, fh, indent=2)

    med = {g: round(stats_blob["groups"][g]["median"], 3) for g in GROUP_ORDER}
    print(f"score_by_group: n={dict((g, stats_blob['groups'][g]['n']) for g in GROUP_ORDER)}")
    print(f"  median rank_pctile: {med}")
    print(f"  Kruskal-Wallis H={H:.1f}, p={p_kw:.3e}")


if __name__ == "__main__":
    main()
