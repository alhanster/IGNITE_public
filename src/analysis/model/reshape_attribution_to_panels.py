#!/usr/bin/env python
"""compute_nb_significance.py must run next: it overwrites panelA/C CSVs with NB-corrected p-values, not Wilcoxon."""
import os
import json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
SRC = os.path.join(ROOT, "outputs", "model", "functional_genomics_contribution",
                   "attribution_results.json")
OUTD = os.path.join(ROOT, "figure_data")
os.makedirs(OUTD, exist_ok=True)


def stars(p):
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2
            else "*" if p < 5e-2 else "ns")


def main():
    d = json.load(open(SRC))
    ca = d["coalition_auc"]
    n_auc = d.get("n_auc_samples", 25)

    def rung(coal_key, model, label, kind):
        m = ca[coal_key]
        return {"model": model, "label": label, "auc": m["mean"], "sd": m["sd"],
                "sem": m["sd"] / (n_auc ** 0.5), "kind": kind}
    panelA = pd.DataFrame([
        # Only arg 3 (label) is editable: 'model' is a key used by _LADDER_KEY and write_submission_readme.py.
        rung("genetic", "genetic only", "Genetics\nOnly", "baseline"),
        rung("genetic+observational", "+ observational", "Genetics\n+ Observational", "single"),
        rung("genetic+perturbational", "+ perturbational", "Genetics\n+ Perturbational", "single"),
        rung("genetic+observational+perturbational", "full", "Full\nModel", "full"),
    ])
    panelA.to_csv(os.path.join(OUTD, "panelA_auc_ladder.csv"), index=False)

    pd_ = d["paired_deltas"]
    panelA_sig = pd.DataFrame([
        {"model_a": "genetic", "model_b": "genetic+observational",
         "delta": pd_["observational_over_genetic"]["delta"],
         "wilcoxon_p": pd_["observational_over_genetic"]["wilcoxon_p"],
         "stars": stars(pd_["observational_over_genetic"]["wilcoxon_p"])},
        {"model_a": "genetic", "model_b": "genetic+perturbational",
         "delta": pd_["perturbational_over_genetic"]["delta"],
         "wilcoxon_p": pd_["perturbational_over_genetic"]["wilcoxon_p"],
         "stars": stars(pd_["perturbational_over_genetic"]["wilcoxon_p"])},
        {"model_a": "genetic", "model_b": "full",
         "delta": pd_["full_vs_genetic"]["delta"],
         "wilcoxon_p": pd_["full_vs_genetic"]["wilcoxon_p"],
         "stars": stars(pd_["full_vs_genetic"]["wilcoxon_p"])},
    ])
    panelA_sig.to_csv(os.path.join(OUTD, "panelA_significance.csv"), index=False)

    # panelB_shapley_shares.csv
    sh_pct, sh_auc = d["shapley_share_pct"], d["shapley_above_chance"]
    panelB = pd.DataFrame([{"evidence": g, "share_pct": sh_pct[g], "shapley_auc": sh_auc[g]}
                           for g in ("genetic", "perturbational", "observational")])
    panelB.to_csv(os.path.join(OUTD, "panelB_shapley_shares.csv"), index=False)

    # panelC_marginal.csv
    panelC = pd.DataFrame([
        {"comparison": "genetic+perturbational\nvs genetic", "evidence": "perturbational",
         "delta": pd_["perturbational_over_genetic"]["delta"], "sem": pd_["perturbational_over_genetic"]["sem"],
         "wilcoxon_p": pd_["perturbational_over_genetic"]["wilcoxon_p"]},
        {"comparison": "genetic+obs\nvs genetic", "evidence": "observational",
         "delta": pd_["observational_over_genetic"]["delta"], "sem": pd_["observational_over_genetic"]["sem"],
         "wilcoxon_p": pd_["observational_over_genetic"]["wilcoxon_p"]},
        {"comparison": "full\nvs genetic", "evidence": "full",
         "delta": pd_["full_vs_genetic"]["delta"], "sem": pd_["full_vs_genetic"]["sem"],
         "wilcoxon_p": pd_["full_vs_genetic"]["wilcoxon_p"]},
    ])
    panelC.to_csv(os.path.join(OUTD, "panelC_marginal.csv"), index=False)

    # panels_ABC_meta.json
    meta = {
        "full_above_chance": d["full_above_chance"],
        "shapley_total": sum(sh_auc.values()),
        "coalition_auc": {
            "chance": ca["chance"]["mean"],
            "genetic": ca["genetic"]["mean"],
            "perturbational": ca["perturbational"]["mean"],
            "observational": ca["observational"]["mean"],
            "genetic+perturbational": ca["genetic+perturbational"]["mean"],
            "genetic+observational": ca["genetic+observational"]["mean"],
            "observational+perturbational": ca["observational+perturbational"]["mean"],
            "genetic+observational+perturbational": ca["genetic+observational+perturbational"]["mean"],
        },
    }
    json.dump(meta, open(os.path.join(OUTD, "panels_ABC_meta.json"), "w"), indent=2)

    print("wrote panelA_auc_ladder.csv, panelA_significance.csv, panelB_shapley_shares.csv,")
    print("      panelC_marginal.csv, panels_ABC_meta.json to", OUTD)


if __name__ == "__main__":
    main()
