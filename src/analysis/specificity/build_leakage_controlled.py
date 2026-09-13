"""Backing data for Supplementary Fig S3: leakage-controlled comparison of GPS, missense-z, and the PU model against held-out trial genes.

Task A ranks every scorable gene and evaluates how well each method separates the 727 held-out in-trial genes from the 18,435 non-target genes, excluding the 340 approved targets used as PU-model training positives. The leakage-controlled arm additionally excludes the 453 genes GPS was trained on, leaving 535 positives untrained on by either method.

Inputs (read-only):
  data/perturbseq/pu/pu_model_matrix.parquet             gene, pu_role, mis.z_score
  outputs/model/full_model_pu_scores.csv                       pu_score (written by pu_target_model.py)
  data/specificity/gps_per_gene.csv                      gene -> gps_max_overall
  data/specificity/gps_
"""
import os, sys, json, shutil
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model"))
# delong_test is imported from compute_delong_specificity.py (Sun-Xu midrank) to avoid a third DeLong implementation beside build_heldout_trial_auc.py's naive one.
from compute_delong_specificity import delong_test
from _stats import holm
import _platform
# pum's feature list and count are imported from pu_target_model.py so they exactly match the features it trains on.
import pu_target_model as pum
import pyarrow.parquet as pq

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
OUTD = os.path.join(ROOT, "figure_data")

MATRIX   = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
PU_PATH  = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
GPS_PATH = os.path.join(ROOT, "data", "specificity", "gps_per_gene.csv")
LAB_PATH = os.path.join(ROOT, "data", "specificity", "gps_drug_labeled_genes.csv")

METHODS = [("PU full model", "pu_score",        "drop"),
           ("GPS overall",   "gps_max_overall", "fill0"),
           ("missense-z",    "mis.z_score",     "drop")]


def stars(p):
    """Significance thresholds: * p<.05, ** p<.01, *** p<.001, **** p<1e-4, computed in this stage rather than the renderer. See REPRODUCIBILITY.md."""
    return ("****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2
            else "*" if p < 5e-2 else "ns")


def main():
    for p in (MATRIX, PU_PATH, GPS_PATH, LAB_PATH):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p}\nRun `make tables` first -- full_model_pu_scores.csv is written by "
                "pu_target_model.py into gitignored outputs/.")

    m   = pd.read_parquet(MATRIX)[["gene", "pu_role", "mis.z_score"]]
    pu  = pd.read_csv(PU_PATH)[["gene", "pu_score"]]
    # gene_id is dropped; each method is scored on its own subframe, so this cannot change any AUC, verified against the committed six-row output file.
    gps = pd.read_csv(GPS_PATH)[["gene", "gps_max_overall"]]
    labeled = set(pd.read_csv(LAB_PATH).iloc[:, 0])

    assert m.gene.is_unique and pu.gene.is_unique, "gene keys must be unique to join on"
    d = (m.merge(pu, on="gene", how="left")
           .merge(gps, on="gene", how="left"))
    assert len(d) == len(m), "joins must not duplicate rows"
    assert d.pu_score.notna().all(), "every gene in the matrix needs a PU score"

    # Task A evaluation set: held-out trial genes vs non-targets; the 340 training positives are excluded, not scored as negatives.
    task = d[d.pu_role != "P"].copy()
    task["y"] = (task.pu_role == "trial_heldout").astype(int)
    ctl = task[~task.gene.isin(labeled)].copy()

    n_lab_pos = int(task.y.sum() - ctl.y.sum())
    print(f"GPS training list: {len(labeled)} genes, {n_lab_pos} of them held-out trial targets")
    print(f"standard: {len(task)} genes / {int(task.y.sum())} positives")
    print(f"leakage-controlled: {len(ctl)} genes / {int(ctl.y.sum())} positives")

    def scored(sub, col, policy):
        """Apply the method's own missing-value policy and return (scores, labels)."""
        s = sub.copy()
        if policy == "fill0":
            s[col] = s[col].fillna(0.0)
        else:
            s = s.dropna(subset=[col])
        return s[col].values.astype(float), s.y.values

    rows = []
    for arm, sub, note in [("standard_taskA", task, "all genes"),
                           ("leakage_controlled", ctl, "GPS drug-labeled genes removed")]:
        for name, col, policy in METHODS:
            sc, y = scored(sub, col, policy)
            rows.append({"method": name, "evaluation": arm,
                         "auc": round(float(roc_auc_score(y, sc)), 3),
                         "n_genes": len(y), "n_pos": int(y.sum()),
                         "missing_policy": policy, "note": note})

    out = pd.DataFrame(rows)
    os.makedirs(OUTD, exist_ok=True)
    out.to_csv(os.path.join(OUTD, "leakage_controlled_comparison.csv"), index=False)
    print(out.to_string(index=False))

    # Copied byte-for-byte via shutil.copyfile, not re-serialized through pandas. See REPRODUCIBILITY.md.
    assert len(labeled) == 453, f"the GPS training list is 453 genes, this one has {len(labeled)}"
    shutil.copyfile(LAB_PATH, os.path.join(OUTD, "gps_drug_labeled_genes.csv"))
    print(f"copied {len(labeled)}-gene GPS training list -> figure_data/gps_drug_labeled_genes.csv")

    # Paired DeLong significance test, both arms. See REPRODUCIBILITY.md.
    def paired_delong(sub):
        pair = sub.dropna(subset=["pu_score"]).copy()
        pair["gps_filled"] = pair.gps_max_overall.fillna(0.0)
        y_pair = pair.y.values
        assert len(pair) == len(sub), "pu_score is never NaN, so the paired set is the whole arm"
        aucs, z, p = delong_test(pair.pu_score.values.astype(float),
                                 pair.gps_filled.values.astype(float), y_pair)
        return {"n_genes": int(len(y_pair)), "n_pos": int(y_pair.sum()),
                "auc_pu": float(aucs[0]), "auc_gps": float(aucs[1]),
                "delta": float(aucs[0] - aucs[1]),
                "z": float(z), "p": float(p)}

    arms = {"standard_taskA": paired_delong(task),
            "leakage_controlled": paired_delong(ctl)}

    # Holm correction across both arms as one family. See REPRODUCIBILITY.md.
    for a_, p_holm in zip(arms.values(), holm([a_["p"] for a_ in arms.values()])):
        a_["p_holm"] = p_holm
        a_["stars"] = stars(p_holm)

    # Asserts the published controlled-arm values (+0.053, z=3.62, p=2.9e-4) so the added arm and dropped join cannot silently change them.
    # Re-pinned for the 24-feature model (gene_burden_score added to the genetic block): auc_pu,
    # z and p moved in the 5th decimal, auc_gps is bit-identical since GPS does not depend on the fit.
    lc = arms["leakage_controlled"]
    assert lc["n_genes"] == 18878 and lc["n_pos"] == 535, \
        f"the controlled arm is 18,878 genes / 535 positives, got {lc['n_genes']} / {lc['n_pos']}"
    # The pins are exact only on the reference platform (Apple Silicon macOS). Elsewhere XGBoost
    # and BLAS floating point drift in the 4th decimal, so a moved value is reported, not fatal.
    # See REPRODUCIBILITY.md.
    for k, want in (("auc_pu", 0.6450824144890129), ("auc_gps", 0.5924324183867029),
                    ("z", 3.6223603841632457), ("p", 0.0002919270292188363)):
        _platform.check(abs(lc[k] - want) < 1e-12,
                        f"controlled-arm {k} moved: {lc[k]!r} vs the published {want!r}")

    # Feature count for labels and Table 11 is derived from the matrix's column schema at read time, not hardcoded, to prevent drift.
    n_feat_full = len([c for c in pq.read_schema(MATRIX).names if c not in pum.META])

    delong = {"comparison": "PU full model vs GPS overall, paired DeLong, per evaluation arm",
              "n_feat_full": n_feat_full,
              "arms": arms,
              "note": ("Paired DeLong on the shared gene set, one entry per arm. Closes gap #1 "
                       "of final_plots/supplementary/README_gps_leakage_controlled.md, which "
                       "flagged that the headline contrast had no significance test. The "
                       "standard arm was added 2026-08-31 so the figure could mark it too.")}
    with open(os.path.join(OUTD, "leakage_controlled_delong.json"), "w") as fh:
        json.dump(delong, fh, indent=1)
    for arm, a in arms.items():
        print(f"\npaired DeLong ({arm}): PU {a['auc_pu']:.4f} vs GPS {a['auc_gps']:.4f}, "
              f"delta = {a['delta']:+.4f}, z = {a['z']:.2f}, p = {a['p']:.3g}  [{a['stars']}]")


if __name__ == "__main__":
    main()
