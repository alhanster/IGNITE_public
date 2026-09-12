#!/usr/bin/env python3
"""Generates panelD_recovery.csv: immune-indication trial recovery at ranking depth,
genetics-only versus the full model, ranked to depth 500.

Membership uses the immune-only rule (furthest immune-indication stage >= Phase 1),
matching build_immune_stages_data.py, so this figure and the trial-validation figure
agree by construction. See REPRODUCIBILITY.md, Other significance tests, for the
ranking provenance of each bar.

Outputs:
  figure_data/panelD_recovery.csv (k, depth, genetics, genetics_perturbseq, gain)
  outputs/model/genetics_only_pu_scores.csv (committed genetics-only ranking)
"""
import os
import numpy as np
import pandas as pd
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"

OUTD      = os.path.join(ROOT, "figure_data")
MODEL_DIR = os.path.join(ROOT, "outputs", "model")
PU_SCORES = os.path.join(MODEL_DIR, "full_model_pu_scores.csv")
OUT_GEN   = os.path.join(MODEL_DIR, "genetics_only_pu_scores.csv")
EVIDENCE  = os.path.join(ROOT, "data", "data_drug", "all_drugs_approved_and_in_trial_evidence_ot.csv")
MATRIX    = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
OUT       = os.path.join(OUTD, "panelD_recovery.csv")

DEPTHS  = (50, 100, 200, 500)                               # uncapped (see docstring)
GENETIC = ["lof.oe_ci.upper", "mis.z_score", "IEI", "gwas_score", "gene_burden_score"]
# Immune-stage ordinal matches build_immune_stages_data.py, required for cross-script comparisons.
STAGE_ORD = {"Approved": 4, "Phase 3": 3, "Phase 2": 2, "Phase 1": 1, "Unknown": 0.5}

# Reuses the production scorer (SEED, T_BAG, pu_bag, META) verbatim so scores match the production run.
_spec = importlib.util.spec_from_file_location(
    "pum", os.path.join(ROOT, "src", "analysis", "pu_target_model.py"))
pum = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pum)


def immune_num_by_gene(ev):
    """Per-gene furthest immune-indication stage (>=1 means an immune trial).

    gene_target splits on ';' (e.g. ruxolitinib -> JAK1;JAK2;JAK3;TYK2) so each gene is
    credited with its drug's stage; otherwise ~2,617/7,030 multi-gene rows file under a
    composite key and undercount recovery.
    """
    ev = ev.copy()
    ev["stg_num"] = ev["drug_stage_for_immune_indication"].map(STAGE_ORD)
    ev = ev.assign(gene_target=ev["gene_target"].str.split(";")).explode("gene_target")
    return ev.groupby("gene_target")["stg_num"].max()


def recovery_counts(gene_order, imm_num, depths=DEPTHS):
    """Reports, for the top-k genes in the given ranked order, the count with an immune trial at stage >= Phase 1."""
    imm = gene_order.map(imm_num).fillna(-1).values >= 1
    return {k: int(imm[:k].sum()) for k in depths}


def main():
    ev = pd.read_csv(EVIDENCE)
    imm_num = immune_num_by_gene(ev)

    # Full-model bar: production scores over the non-approved pool, corresponds to Figure 5b.
    g = pd.read_csv(PU_SCORES)
    nonappr = g[g.pu_role != "P"].sort_values("pu_score", ascending=False).reset_index(drop=True)
    full = recovery_counts(nonappr.gene, imm_num)

    # Genetics-only bar: re-scores the same pool using the production configuration.
    m = pd.read_parquet(MATRIX)
    role = m.pu_role.values
    P_idx = np.where(role == "P")[0]
    U_idx = np.where(role == "unlabeled")[0]
    X = m[GENETIC].apply(pd.to_numeric, errors="coerce").values.astype(float)
    sc, _, _ = pum.pu_bag(X, P_idx, U_idx, GENETIC, T=pum.T_BAG, seed=pum.SEED)
    m_na = m[m.pu_role != "P"].copy()
    m_na["sc"] = sc[role != "P"]
    m_na = m_na.sort_values("sc", ascending=False).reset_index(drop=True)
    gen = recovery_counts(m_na.gene, imm_num)

    # Persists the genetics-only ranking; see REPRODUCIBILITY.md, Other significance tests.
    cols = ["gene", "pu_role", "furthest_stage"]
    ranked = m_na[cols].assign(genetics_score=m_na["sc"],
                               rank=np.arange(1, len(m_na) + 1))
    appr = m[m.pu_role == "P"][cols].assign(genetics_score=sc[role == "P"], rank=pd.NA)
    gen_scores = pd.concat([ranked, appr], ignore_index=True)[
        ["gene", "genetics_score", "pu_role", "furthest_stage", "rank"]]
    # Rank uses nullable Int64 so approved genes (no rank) store as an empty field, not float64's '1.0'.
    gen_scores["rank"] = gen_scores["rank"].astype("Int64")
    os.makedirs(MODEL_DIR, exist_ok=True)
    gen_scores.to_csv(OUT_GEN, index=False)

    out = pd.DataFrame([{"k": k, "depth": f"top-{k}",
                         "genetics": gen[k], "genetics_perturbseq": full[k],
                         "gain": full[k] - gen[k]} for k in DEPTHS])
    os.makedirs(OUTD, exist_ok=True)
    out.to_csv(OUT, index=False)
    print("wrote", OUT)
    print("wrote", OUT_GEN, f"({len(gen_scores)} rows, {len(ranked)} ranked)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
