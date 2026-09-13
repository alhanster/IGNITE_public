#!/usr/bin/env python3
"""build_drug_status_table.py: per-gene immune-indication drug annotations, giving the label, furthest trial phase, the drugs and indications behind it, and trial or approval years.

Step of `make tables`. Reads only committed inputs and joins them: no fitting, no seed, no network. Independent of pu_target_model.py.

Inputs:
  data/data_drug/all_drugs_approved_and_in_trial_evidence_ot.csv   per-drug evidence, the authoritative grain: 7,030 rows, 1,611 drugs, gene_target ';'-joined for multi-target drugs
  data/data_drug/all_drugs_approved_and_in_trial_by_gene_ot.csv    per-gene rollup, used for cross-checks
  data/prospective/drug_trial_dates.json       per-drug CT.gov first-trial years
  data/prospective/positives_dated.csv         approval year for the 340 approved genes
  data/perturbseq/pu/pu_model_matrix.parquet   the 19,502-gene modelled universe

Output:
  figure_data/drug_status_annotations.csv      1,089 rows, one per gene with any annotation

See REPRODUCIBILITY.md, Drug status table, for why list columns come from the evidence file rather than the rollup and why latest_trial_year is omitted.

Run: PYTHONPATH=src .venv/bin/python src/analysis/tables/build_drug_status_table.py
"""
import json
import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"

EVIDENCE = os.path.join(ROOT, "data", "data_drug",
                        "all_drugs_approved_and_in_trial_evidence_ot.csv")
BY_GENE = os.path.join(ROOT, "data", "data_drug",
                       "all_drugs_approved_and_in_trial_by_gene_ot.csv")
DATES = os.path.join(ROOT, "data", "prospective", "drug_trial_dates.json")
APPROVALS = os.path.join(ROOT, "data", "prospective", "positives_dated.csv")
MATRIX = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
OUT = os.path.join(ROOT, "figure_data", "drug_status_annotations.csv")

# Open Targets release version: see PROVENANCE.md.
OT_RELEASE = "26.06"

# Most-advanced-first; Unknown sorts last, keeping its 15 genes in the unlabeled pool, not the trial set.
STAGE_ORDER = ["Approved", "Phase 3", "Phase 2", "Phase 1", "Unknown"]


def drug_status(furthest_stage):
    """Matches build_score_by_group.py: pu_role=='P' equals the 340 Approved-stage genes, asserted in code."""
    if furthest_stage == "Approved":
        return "approved"
    return "in-trial" if str(furthest_stage).startswith("Phase") else "non-target"


def main():
    ev = pd.read_csv(EVIDENCE)

    # Explodes ';'-joined multi-target drugs to one row per (drug, gene); MT-3724 yields 79 rows (PROVENANCE.md).
    ev["gene"] = ev["gene_target"].astype(str).str.split(";")
    ev = ev.explode("gene")
    ev["gene"] = ev["gene"].str.strip()
    ev = ev[ev["gene"] != ""]

    stage_rank = {s: i for i, s in enumerate(STAGE_ORDER)}
    unknown = set(ev["drug_stage_for_immune_indication"]) - set(stage_rank)
    assert not unknown, f"unhandled trial stage(s): {sorted(unknown)}"
    ev["_rank"] = ev["drug_stage_for_immune_indication"].map(stage_rank)

    # Reports min(any_year) per gene, the first trial in any condition; immune_year is null more often, hence the name.
    with open(DATES) as fh:
        dates = json.load(fh)
    year_of = {k: v.get("any_year") for k, v in dates.items() if isinstance(v, dict)}
    missing_dates = set(ev["chembl_id"]) - set(year_of)
    assert not missing_dates, f"{len(missing_dates)} drugs absent from the CT.gov date cache"
    ev["_year"] = ev["chembl_id"].map(year_of)

    g = ev.groupby("gene")
    out = pd.DataFrame({
        "n_drugs": g["chembl_id"].nunique(),
        "drug_names": g["drug"].apply(lambda s: "; ".join(sorted(set(s)))),
        "indications": g["immune_indication"].apply(lambda s: "; ".join(sorted(set(s.dropna())))),
        "furthest_stage": g["_rank"].min().map(dict(enumerate(STAGE_ORDER))),
        "first_trial_year": g["_year"].min(),
    }).reset_index()
    out["drug_status"] = out["furthest_stage"].map(drug_status)

    # Expanded means no single-target drug reaches the gene; it enters via a multi-target drug only.
    single = set(ev.loc[~ev["gene_target"].astype(str).str.contains(";"), "gene"])
    out["multi_target_expanded"] = ~out["gene"].isin(single)

    appr = pd.read_csv(APPROVALS)
    assert appr["gene"].is_unique
    out = out.merge(appr[["gene", "approval_year"]], on="gene", how="left", validate="one_to_one")

    matrix_genes = set(pd.read_parquet(MATRIX, columns=["gene"])["gene"])
    out["in_feature_matrix"] = out["gene"].isin(matrix_genes)

    out = out[["gene", "drug_status", "furthest_stage", "n_drugs", "drug_names", "indications",
               "first_trial_year", "approval_year", "multi_target_expanded",
               "in_feature_matrix"]].sort_values("gene").reset_index(drop=True)

    # Cross-checks against the rollup and model matrix: see REPRODUCIBILITY.md, Build-time assertions inside analysis scripts.
    ref = pd.read_csv(BY_GENE).rename(columns={"gene_target": "gene"})
    assert set(ref["gene"]) == set(out["gene"]), "gene set disagrees with the by-gene rollup"
    chk = out.merge(ref[["gene", "n_drugs", "furthest_stage"]], on="gene", suffixes=("", "_ref"))
    bad_n = chk.loc[chk["n_drugs"] != chk["n_drugs_ref"], "gene"]
    assert bad_n.empty, f"n_drugs disagrees with the rollup for {len(bad_n)} genes"
    bad_s = chk.loc[chk["furthest_stage"] != chk["furthest_stage_ref"], "gene"]
    assert bad_s.empty, f"furthest_stage disagrees with the rollup for {len(bad_s)} genes"

    # Approved set is exactly pu_role=='P', kept equivalent to assign_group's pu_role-first branch.
    roles = pd.read_parquet(MATRIX, columns=["gene", "pu_role"])
    p_genes = set(roles.loc[roles["pu_role"] == "P", "gene"])
    assert set(out.loc[out["drug_status"] == "approved", "gene"]) == p_genes, \
        "approved set is not exactly pu_role == 'P'"

    # Genes in trial but absent from the matrix are the mitochondrial NADH subunits, named not just counted.
    off = set(out.loc[(out["drug_status"] == "in-trial") & ~out["in_feature_matrix"], "gene"])
    assert off == {"MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"}, \
        f"unexpected in-trial genes outside the feature matrix: {sorted(off)}"
    n_trial_in_matrix = int(((out["drug_status"] == "in-trial") & out["in_feature_matrix"]).sum())
    assert n_trial_in_matrix == 727, f"expected 727 trial_heldout genes, got {n_trial_in_matrix}"

    # Renamed only here, after the rollup cross-check above compares it against that file's own
    # furthest_stage. Published name says which indications the stage is for; see build_ranked_atlas_table.py.
    out = out.rename(columns={"furthest_stage": "furthest_imm_stage"})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}")
    print(f"  {len(out)} genes | Open Targets {OT_RELEASE} | "
          f"status: {out['drug_status'].value_counts().to_dict()}")
    print(f"  approval_year populated for {int(out['approval_year'].notna().sum())} of "
          f"{int((out['drug_status'] == 'approved').sum())} approved; "
          f"first_trial_year null for {int(out['first_trial_year'].isna().sum())}; "
          f"multi_target_expanded {int(out['multi_target_expanded'].sum())}")


if __name__ == "__main__":
    main()
