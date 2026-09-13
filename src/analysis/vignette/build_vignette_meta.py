#!/usr/bin/env python3
"""build_vignette_meta.py

Writes figure_data/vignette_meta.json, the panel a and c values used by fig_stat4_vignette.R.

Inputs: outputs/model/full_model_pu_scores.csv (stat4_score, stat4_rank, approved_median,
n_genes, stat4_top_pct); figure_data/vignette_knn_stats.json (knn_p95, knn_median_unlabeled);
figure_data/vignette_panelC_gwas_meta.json (ot_diseases, approved_drugs). gwas_total is
carried from the seed file; see PROVENANCE.md.

Run after pu_target_model regenerates full_model_pu_scores.csv:
    python src/analysis/vignette/build_vignette_meta.py
"""
import os
import json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SCORES = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
# SEED_META is the read-only committed seed; derived fields are merged onto a copy, not written back in place.
SEED_META = os.path.join(ROOT, "data", "case_vignette", "vignette_meta.json")
OUT_META = os.path.join(ROOT, "figure_data", "vignette_meta.json")
KNN_STATS = os.path.join(ROOT, "figure_data", "vignette_knn_stats.json")
# ot_diseases and approved_drugs are written by the opt-in fetch_vignette_panelC_gwas.py script.
PANELC_META = os.path.join(ROOT, "figure_data", "vignette_panelC_gwas_meta.json")
FOCAL = "STAT4"
MODEL_FIELDS = ("stat4_score", "stat4_rank", "approved_median")
# Local name differs from the PANELC_META key: named for content (OT associations for STAT4), not use.
PANELC_FIELDS = {"ot_diseases": "n_diseases_total", "approved_drugs": "approved_drugs"}
# See PROVENANCE.md, STAT4 vignette — `final_plots/figure_stat4_vignette.png`.
PRESERVED = ("gwas_total",)
# Key order is pinned explicitly: json.dump keeps insertion order and verify-tables byte-compares it.
KEY_ORDER = ("stat4_score", "approved_median", "stat4_rank", "knn_p95", "gwas_total",
             "ot_diseases", "approved_drugs", "knn_median_unlabeled", "n_genes",
             "stat4_top_pct")


def main():
    g = pd.read_csv(SCORES)
    if FOCAL not in set(g.gene):
        raise SystemExit("%s not found in %s" % (FOCAL, SCORES))

    focal = g[g.gene == FOCAL].iloc[0]
    approved_median = float(g[g.pu_role == "P"]["pu_score"].median())
    ordered = g.sort_values("pu_score", ascending=False).reset_index(drop=True)
    rank = int(ordered.index[ordered.gene == FOCAL][0]) + 1

    if not os.path.exists(SEED_META):
        raise SystemExit(
            "missing %s -- it carries gwas_total, which nothing in this repo "
            "recomputes" % SEED_META)
    meta = json.load(open(SEED_META))

    # See REPRODUCIBILITY.md, Build-time assertions inside analysis scripts.
    for pth, why in ((KNN_STATS, "knn_p95, panel d's reference line"),
                     (PANELC_META, "ot_diseases / approved_drugs, no longer drawn but "
                                   "still byte-verified in vignette_meta.json")):
        if not os.path.exists(pth):
            raise SystemExit("missing %s -- it carries %s" % (pth, why))
    knn = json.load(open(KNN_STATS))
    meta["knn_p95"] = knn["knn_p95"]
    meta["knn_median_unlabeled"] = knn["knn_median_unlabeled"]

    panelc = json.load(open(PANELC_META))
    for k, src in PANELC_FIELDS.items():
        meta[k] = panelc[src]
    meta["stat4_score"] = round(float(focal.pu_score), 6)
    meta["approved_median"] = round(approved_median, 6)
    meta["stat4_rank"] = rank
    # n_genes and stat4_top_pct are the rank's denominator and panel a's top-percent value for STAT4.
    meta["n_genes"] = int(len(g))
    meta["stat4_top_pct"] = 100.0 * rank / len(g)

    missing = [k for k in KEY_ORDER if k not in meta]
    extra = [k for k in meta if k not in KEY_ORDER]
    if missing or extra:
        raise SystemExit("KEY_ORDER does not describe what was built: missing %s, extra %s"
                         % (missing, extra))
    os.makedirs(os.path.dirname(OUT_META), exist_ok=True)
    json.dump({k: meta[k] for k in KEY_ORDER}, open(OUT_META, "w"), indent=1)

    print("wrote %s (computed from %s):"
          % (os.path.relpath(OUT_META, ROOT), os.path.relpath(SCORES, ROOT)))
    for k in MODEL_FIELDS:
        print("  %-16s %s" % (k, meta[k]))
    print("  (%s rank %d / %d = top %.3f%%)" % (FOCAL, rank, len(g), 100 * rank / len(g)))
    print("  %-16s %s  (from vignette_knn_stats.json)" % ("knn_p95", meta["knn_p95"]))
    print("  panel c counts from Open Targets %s via fetch_vignette_panelC_gwas.py: %s"
          % (panelc["ot_release"], {k: meta[k] for k in PANELC_FIELDS}))
    print("carried over from the seed (no producer anywhere):",
          {k: meta[k] for k in PRESERVED})


if __name__ == "__main__":
    main()
