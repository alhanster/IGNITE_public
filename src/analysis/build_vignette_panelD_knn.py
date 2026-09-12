"""STAT4 vignette panel d: nearest approved-target neighbors in knockdown-signature space.

Scores each approved immune-drug target by the cosine similarity of its transcriptome-wide
knockdown signature to STAT4's, in the stimulated 48-hour condition, and reports the top six
approved targets against the unlabeled distribution. This is an annotation for a case
vignette, not a screen: held-out trial-versus-unlabeled discrimination from this cosine is
near chance (AUC 0.47-0.53).

Inputs: data/perturbseq/knn/signatures_Stim48hr.npy (knockdown signatures, not committed, rebuilt by
`make knn-signatures`), matching gene symbol files, the HGNC symbol table, and
outputs/model/full_model_pu_scores.csv for the anchor set (pu_role == "P").

Outputs: figure_data/vignette_panelD_knn.csv (six neighbors: gene, cosine, mechanism,
pathway, approved drugs) and vignette_knn_stats.json (knn_p95).

See REPRODUCIBILITY.md, knn-signatures, for the signature-matrix skip
behavior and verification values.
"""
import os, sys, json
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src", "perturbseq_features"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gene_name_utils as gnu
from _skip_ledger import record_skip  # noqa: E402

FIGDATA = os.path.join(ROOT, "figure_data")
KNN_DATA = os.path.join(ROOT, "data", "perturbseq", "knn")
HGNC = os.path.join(ROOT, "data", "perturbseq", "sources", "hgnc_symbol_subset.txt")
SCORES = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
# Approved-drug lists are committed static inputs; unlike SCORES, they impose no ordering requirement.
DRUGS_BY_GENE = os.path.join(ROOT, "data", "data_drug",
                             "all_drugs_approved_and_in_trial_by_gene_ot.csv")
DRUG_EVIDENCE = os.path.join(ROOT, "data", "data_drug",
                             "all_drugs_approved_and_in_trial_evidence_ot.csv")

FOCAL, COND, TOPN = "STAT4", "Stim48hr", 6
MAX_DRUGS_SHOWN = 1          # then "+ N more"

# Hand-curated columns, not derived from other data; see PROVENANCE.md.
ANNOTATION = {
    "CD3E":   ("CD3 complex",         "cd3"),
    "LCK":    ("TCR proximal kinase", "tcrk"),
    "CD3G":   ("CD3 complex",         "cd3"),
    "CD2":    ("co-stim / adhesion",  "costim"),
    "PPP3CA": ("calcineurin",         "calci"),
    "PRKCQ":  ("PKCθ (TCR→NF-κB)", "tcrk"),
}


def _pretty(name):
    """Hyphen parts with a digit, or very short, stay uncapitalized (e.g. MUROMONAB-CD3 -> Muromonab-CD3)."""
    def part(p):
        return p if (any(c.isdigit() for c in p) or len(p) <= 3) else p.title()
    return " ".join("-".join(part(p) for p in tok.split("-")) for tok in name.split())


def approved_drug_columns():
    """Maps gene to (all approved drugs joined by "; ", the truncated panel label), using furthest_stage_drugs

    (approved only). See REPRODUCIBILITY.md, Vignette metadata.
    """
    bg = pd.read_csv(DRUGS_BY_GENE)
    ev = pd.read_csv(DRUG_EVIDENCE)
    immune = {g: set(d.drug.astype(str).str.upper()) for g, d in ev.groupby("gene_target")}
    cols = {}
    for _, r in bg.iterrows():
        approved = [d for d in str(r.furthest_stage_drugs).split(";") if d and d != "nan"]
        if not approved:
            continue
        imm = immune.get(r.gene_target, set())
        ordered = sorted(approved, key=lambda d: (d not in imm, d))
        # One ordering and _pretty() serve both sheet and panel, preventing drift between them.
        pretty = [_pretty(d) for d in ordered]
        shown = " \u00b7 ".join(pretty[:MAX_DRUGS_SHOWN])
        extra = len(pretty) - MAX_DRUGS_SHOWN
        cols[r.gene_target] = ("; ".join(pretty),
                               shown + (f"  + {extra} more" if extra > 0 else ""))
    return cols


def anchor_set():
    """The 340 approved immune targets, HGNC-harmonized (no-op today, kept for future symbol changes)."""
    scores = pd.read_csv(SCORES)
    approved = scores.loc[scores["pu_role"] == "P", "gene"].tolist()
    hgnc = gnu.load_hgnc(path=HGNC)
    out = gnu.harmonize(pd.DataFrame({"gene": approved}), hgnc)
    ref = set((out[0] if isinstance(out, tuple) else out)["gene"])
    prev2app, approved_sym = hgnc["prev2app"], hgnc["approved"]
    canon = lambda x: prev2app[x] if (x not in approved_sym and x in prev2app) else x
    return ref, canon


def main():
    # See REPRODUCIBILITY.md, knn-signatures, for the skip behavior on a fresh clone.
    sig = os.path.join(KNN_DATA, f"signatures_{COND}.npy")
    if not os.path.exists(sig):
        print(f"SKIP {os.path.basename(__file__)}: {os.path.relpath(sig, ROOT)} is absent.\n"
              "  It is not committed (294 MB). Run `make knn-signatures` to stream it from the\n"
              "  public source (~3-4 min), then re-run this script to re-derive panel d.\n"
              "  figure_data/vignette_panelD_knn.csv and vignette_knn_stats.json are committed\n"
              "  and left untouched, so every downstream step and `make verify-tables` are fine.")
        record_skip(os.path.basename(__file__),
                    "signature matrix %s absent" % os.path.basename(sig),
                    ["figure_data/vignette_panelD_knn.csv",
                     "figure_data/vignette_knn_stats.json"])
        return

    # Missing full_model_pu_scores.csv or the HGNC table hard-fails; both are always-expected inputs.
    for p in (SCORES, HGNC):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p}\nfull_model_pu_scores.csv comes from `make tables` -- run it first.")

    reference, canon = anchor_set()
    X = np.load(sig).astype(np.float32)
    rg = np.loadtxt(os.path.join(KNN_DATA, f"signatures_{COND}.genes.txt"), dtype=str)
    var_name = np.loadtxt(os.path.join(KNN_DATA, "var_gene_names.txt"), dtype=str)
    print(f"signatures {X.shape} ({COND}); {len(reference)} approved-target anchors")

    # Each gene's column is zeroed before normalising, so its knockdown can't inflate its own score.
    name2col = {g: i for i, g in enumerate(var_name)}
    for r, g in enumerate(rg):
        c = name2col.get(g)
        if c is not None:
            X[r, c] = 0.0
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)

    gi = np.where(rg == FOCAL)[0]
    if not len(gi):
        raise SystemExit(f"{FOCAL} has no knockdown signature in {COND}")
    sims = Xn @ Xn[gi[0]]
    sims[gi[0]] = -np.inf                       # never a neighbor of itself

    is_anchor = np.array([canon(g) in reference for g in rg])
    # STAT4's -inf sentinel stays in the unlabeled pool rather than being special-cased; see
    #     REPRODUCIBILITY.md, Vignette metadata.
    unl = sims[~is_anchor]
    pos_g, pos_v = rg[is_anchor], sims[is_anchor]
    order = np.argsort(pos_v)[::-1][:TOPN]

    drug_cols = approved_drug_columns()
    rows = []
    for rank, i in enumerate(order, start=1):
        g = str(pos_g[i])
        mech, path = ANNOTATION.get(g, ("immune target", "other"))
        full, label = drug_cols.get(g, ("", ""))
        rows.append({"gene": g, "cosine": round(float(pos_v[i]), 3),
                     "mechanism": mech, "approved_drugs": full,
                     "approved_drugs_label": label, "pathway": path,
                     "order": rank})
    out = pd.DataFrame(rows)

    nodrugs = [r["gene"] for r in rows if not r["approved_drugs"]]
    if nodrugs:
        # An empty approved-drug list means the gene is missing from the OT rollup, not drug-free.
        print(f"WARNING: no approved-drug row for {nodrugs} in "
              f"{os.path.basename(DRUGS_BY_GENE)} -- the panel would draw a blank line",
              file=sys.stderr)

    unannotated = [r["gene"] for r in rows if r["gene"] not in ANNOTATION]
    if unannotated:
        # Fails loudly, rather than rendering blank mechanisms, if curated labels miss changed neighbors.
        print(f"WARNING: no curated annotation for {unannotated} -- "
              f"update ANNOTATION in {os.path.basename(__file__)}", file=sys.stderr)

    os.makedirs(FIGDATA, exist_ok=True)
    out.to_csv(os.path.join(FIGDATA, "vignette_panelD_knn.csv"), index=False)
    print(out.to_string(index=False))

    # knn_p95 (panel d's dashed reference line) is read by build_vignette_meta.py; see PROVENANCE.md.
    stats = {"knn_p95": round(float(np.quantile(unl, 0.95)), 3),
             "knn_median_unlabeled": round(float(np.median(unl)), 3),
             "condition": COND, "focal": FOCAL,
             "n_anchor_rows": int(is_anchor.sum()), "n_unlabeled_rows": int((~is_anchor).sum())}
    with open(os.path.join(FIGDATA, "vignette_knn_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
