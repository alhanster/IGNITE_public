#!/usr/bin/env python
"""Stage-1 tables for final_plots/genetic_data.png panels b and c.

Inputs: data/full_gene_list.tsv, data/data_drug/approved_target_genes.txt,
data/data_drug/all_drugs_approved_and_in_trial_by_gene_ot.csv
Outputs: figure_data/gwas_by_drug_status.csv (gene, group, gwas_score; gwas_score>0 only),
figure_data/mis_z_by_drug_status.csv (gene, group, mis_z_score; non-NA),
figure_data/genetic_data_group_n.csv (panel, group, n),
figure_data/gwas_zero_inflation.csv (group, pct_gt0)

See REPRODUCIBILITY.md, Drug status table, for the drug-status
partition rule and the per-panel NA handling.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
DATA = os.path.join(ROOT, "data")
FIGDATA = os.path.join(ROOT, "figure_data")

GENE_LIST    = os.path.join(DATA, "full_gene_list.tsv")
APPROVED_TXT = os.path.join(DATA, "data_drug", "approved_target_genes.txt")
DRUG_CSV     = os.path.join(DATA, "data_drug",
                            "all_drugs_approved_and_in_trial_by_gene_ot.csv")

# Renderer's display order; the fill palette is keyed on it.
TARGET_LEVELS = ["approved", "in_trials", "other"]
TRIAL_STAGES = ["Phase 1", "Phase 2", "Phase 3"]


def drug_status_sets():
    """Approved and trial-only gene sets, matching how the renderer builds them."""
    with open(APPROVED_TXT) as fh:
        approved = {ln.strip() for ln in fh if ln.strip()}

    drug = pd.read_csv(DRUG_CSV)
    drug["gene_target"] = drug["gene_target"].astype(str).str.strip()
    staged = drug.loc[drug["furthest_stage"].isin(TRIAL_STAGES), "gene_target"]
    staged = {g for g in staged if g}
    # Genes marked approved are excluded from the trial group, so each gene is counted once.
    return approved, staged - approved


def assign_group(genes, approved, trial_only):
    """Group assignment follows approved > in_trials > other, matching the renderer's case_when order."""
    return [
        "approved" if g in approved else "in_trials" if g in trial_only else "other"
        for g in genes
    ]


def main():
    approved, trial_only = drug_status_sets()
    fg = pd.read_csv(GENE_LIST, sep="\t")

    # Panel b: immune-GWAS score
    gw = fg[["gene", "gwas_score"]].copy()
    gw["gwas_score"] = pd.to_numeric(gw["gwas_score"], errors="coerce")
    gw = gw[gw["gwas_score"].notna()].copy()
    gw["group"] = assign_group(gw["gene"], approved, trial_only)

    # Counts and the zero-inflation percentage are computed from the unfiltered data frame.
    gw_n = gw.groupby("group").size()
    gw_pct = gw.assign(pos=gw["gwas_score"] > 0).groupby("group")["pos"].mean()

    gw_pos = gw[gw["gwas_score"] > 0][["gene", "group", "gwas_score"]]

    # Panel c: missense-constraint z-score
    # Its own NA filter, so its own counts.
    mz = fg[["gene", "mis.z_score"]].copy()
    mz["mis_z_score"] = pd.to_numeric(mz["mis.z_score"], errors="coerce")
    mz = mz[mz["mis_z_score"].notna()].copy()
    mz["group"] = assign_group(mz["gene"], approved, trial_only)
    mz_n = mz.groupby("group").size()
    mz_out = mz[["gene", "group", "mis_z_score"]]

    # Write output tables.
    os.makedirs(FIGDATA, exist_ok=True)
    gw_pos.to_csv(os.path.join(FIGDATA, "gwas_by_drug_status.csv"), index=False)
    mz_out.to_csv(os.path.join(FIGDATA, "mis_z_by_drug_status.csv"), index=False)

    counts = pd.DataFrame(
        [{"panel": "gwas", "group": g, "n": int(gw_n.get(g, 0))} for g in TARGET_LEVELS]
        + [{"panel": "mis_z", "group": g, "n": int(mz_n.get(g, 0))} for g in TARGET_LEVELS]
    )
    counts.to_csv(os.path.join(FIGDATA, "genetic_data_group_n.csv"), index=False)

    pd.DataFrame(
        [{"group": g, "pct_gt0": float(gw_pct.get(g, 0.0))} for g in TARGET_LEVELS]
    ).to_csv(os.path.join(FIGDATA, "gwas_zero_inflation.csv"), index=False)

    print("genetic_data tables written to figure_data/")
    for g in TARGET_LEVELS:
        print("  %-10s gwas n=%-6d (%.1f%% > 0)   mis.z n=%d"
              % (g, gw_n.get(g, 0), 100 * gw_pct.get(g, 0.0), mz_n.get(g, 0)))


if __name__ == "__main__":
    main()
