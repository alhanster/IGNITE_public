"""Adds approved-drug annotations to each kNN nearest-target table, producing Supplementary Table 15.

Reads outputs/knn/knn_nearest_target_<COND>.csv from build_knn_nearest_target.py and writes figure_data/knn_nearest_target_<COND>.csv (8 columns, one sheet per condition).

See REPRODUCIBILITY.md for the skip behavior.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from _skip_ledger import record_skip  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"

KNN_DIR = os.path.join(ROOT, "outputs", "knn")
FIGDATA = os.path.join(ROOT, "figure_data")
DRUGS = os.path.join(ROOT, "data", "data_drug",
                     "all_drugs_approved_and_in_trial_by_gene_ot.csv")
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
TARGET_COL = "nearest_immune_target"


def build_maps(path):
    """Maps gene to approved drug names and gene to immune indications, both joined with "; " into single strings, not lists."""
    d = pd.read_csv(path)
    ap = d[d["furthest_stage"] == "Approved"]
    dmap, imap = {}, {}
    for _, r in ap.iterrows():
        g = r["gene_target"]
        drugs = str(r["furthest_stage_drugs"]) if pd.notna(r["furthest_stage_drugs"]) else ""
        indic = str(r["immune_indications"]) if pd.notna(r["immune_indications"]) else ""
        dmap[g] = "; ".join(sorted(set(x for x in drugs.split(";") if x)))
        imap[g] = "; ".join(sorted(set(x for x in indic.split(";") if x)))
    return dmap, imap


def annotate(cond, dmap, imap):
    knn = pd.read_csv(os.path.join(KNN_DIR, f"knn_nearest_target_{cond}.csv"))
    tgt = knn[TARGET_COL]

    # Reverse-inserted at one position; final order: nearest_immune_target, drug_name, immune_indications.
    cols = [("immune_indications", tgt.map(lambda g: imap.get(g, ""))),
            ("drug_name", tgt.map(lambda g: dmap.get(g, "")))]
    at = knn.columns.get_loc(TARGET_COL) + 1
    for name, vals in cols:
        if name in knn.columns:
            knn = knn.drop(columns=[name])
        knn.insert(at, name, vals)

    out = os.path.join(FIGDATA, f"knn_nearest_target_{cond}.csv")
    knn.to_csv(out, index=False)
    n = int(tgt.map(lambda g: g in dmap).sum())
    # Missing matches are surfaced, not asserted, since the anchor set comes from pu_role and can drift from the Open Targets rollup.
    print(f"[{cond}] rows={len(knn)}, nearest neighbor has an approved drug: {n}")
    return len(knn)


def main():
    absent = [c for c in CONDS
              if not os.path.isfile(os.path.join(KNN_DIR, f"knn_nearest_target_{c}.csv"))]
    if absent:
        print(f"SKIP {os.path.basename(__file__)}: outputs/knn/ has no table for "
              f"{', '.join(absent)}.\n"
              "  build_knn_nearest_target.py skips when the signature matrices are absent,\n"
              "  which is the normal state of a fresh clone -- run `make knn-signatures`, then\n"
              "  both steps, to re-derive Supplementary Table 15. The three\n"
              "  figure_data/knn_nearest_target_*.csv are committed and left untouched.")
        record_skip(os.path.basename(__file__),
                    "upstream knn tables absent for %s" % ", ".join(absent),
                    ["figure_data/knn_nearest_target_%s.csv" % c for c in CONDS])
        return

    dmap, imap = build_maps(DRUGS)
    print(f"approved-stage genes in the Open Targets rollup: {len(dmap)}")
    os.makedirs(FIGDATA, exist_ok=True)
    total = sum(annotate(c, dmap, imap) for c in CONDS)
    print(f"wrote {len(CONDS)} tables ({total} rows) to figure_data/")
    print("confirm nothing moved: make verify-tables")


if __name__ == "__main__":
    main()
