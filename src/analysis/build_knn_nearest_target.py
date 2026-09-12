"""Nearest approved immune-target knockdown neighbor, per stimulation condition.

For every gene with a significant on-target CD4+ T-cell knockdown, computes the cosine similarity of its
transcriptome-wide signature to the nearest approved immune-drug target anchor:

    score(gene g) = max over anchors t != g of cosine(signature(g), signature(t))

Run for Rest, Stim8hr, and Stim48hr conditions; reported as Supplementary Table 15.

See REPRODUCIBILITY.md, knn-signatures.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"

sys.path.insert(0, os.path.join(ROOT, "src", "perturbseq_features"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gene_name_utils as gnu
from _skip_ledger import record_skip  # noqa: E402

KNN_DATA = os.path.join(ROOT, "data", "perturbseq", "knn")
HGNC = os.path.join(ROOT, "data", "perturbseq", "sources", "hgnc_symbol_subset.txt")
SCORES = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
OUTDIR = os.path.join(ROOT, "outputs", "knn")
CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def signature_paths(cond):
    """(npy, genes.txt) for a condition. Existence is checked once, in main()."""
    return (os.path.join(KNN_DATA, f"signatures_{cond}.npy"),
            os.path.join(KNN_DATA, f"signatures_{cond}.genes.txt"))


def missing_signatures():
    """Repo-relative paths of every signature file absent from data/perturbseq/knn/."""
    return [os.path.relpath(p, ROOT)
            for cond in CONDS for p in signature_paths(cond) if not os.path.isfile(p)]


def anchor_set():
    """340 approved immune targets (pu_role == P), harmonized to HGNC symbols so future symbol changes cannot silently mis-key anchors."""
    scores = pd.read_csv(SCORES)
    approved = scores.loc[scores["pu_role"] == "P", "gene"].tolist()
    hgnc = gnu.load_hgnc(path=HGNC)
    out = gnu.harmonize(pd.DataFrame({"gene": approved}), hgnc)
    ref = set((out[0] if isinstance(out, tuple) else out)["gene"])
    prev2app, approved_sym = hgnc["prev2app"], hgnc["approved"]
    canon = lambda x: prev2app[x] if (x not in approved_sym and x in prev2app) else x
    return ref, canon


def score_condition(cond, reference, canon, var_name):
    npy, genes_txt = signature_paths(cond)
    X = np.load(npy).astype(np.float32)
    rg = np.loadtxt(genes_txt, dtype=str)

    # Zeroes each perturbed gene's own column before normalizing, so its own knockdown effect does not inflate its similarity row; part two is the -inf diagonal below.
    name2col = {g: i for i, g in enumerate(var_name)}
    for r, g in enumerate(rg):
        c = name2col.get(g)
        if c is not None:
            X[r, c] = 0.0
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)

    S = Xn @ Xn.T
    np.fill_diagonal(S, -np.inf)                    # never its own nearest neighbor

    rg_c = np.array([canon(g) for g in rg])
    pos = np.array([g in reference for g in rg_c])
    pix = np.where(pos)[0]
    sc = np.empty(len(rg))
    nn = np.empty(len(rg), dtype=object)
    for i in range(len(rg)):
        cols = pix[pix != i]                        # an anchor is scored against the OTHERS
        j = cols[np.argmax(S[i, cols])]
        sc[i] = S[i, j]
        nn[i] = rg[j]

    # Signature strength, computed on the self-zeroed X, so a gene's own knockdown never inflates its own reported strength.
    strength = (np.abs(X) > 3).sum(1)

    return pd.DataFrame({
        "gene": rg, "kNN_max_score": sc, "nearest_immune_target": nn,
        "is_approved_target": pos, "n_downstream_z3": strength,
        "culture_condition": cond,
    }).sort_values("kNN_max_score", ascending=False).reset_index(drop=True)


def main():
    absent = missing_signatures()
    if absent:
        print(f"SKIP {os.path.basename(__file__)}: {len(absent)} signature file(s) absent, "
              f"first is {absent[0]}.\n"
              "  They are not committed (294 MB each). Run `make knn-signatures` to stream them\n"
              "  from the public source (~3-4 min), then re-run this script and\n"
              "  add_knn_approved_drugs.py to re-derive Supplementary Table 15.\n"
              "  The three figure_data/knn_nearest_target_*.csv are committed and left\n"
              "  untouched, so every downstream step and `make verify-tables` are fine.")
        record_skip(os.path.basename(__file__),
                    "%d signature matrix file(s) absent" % len(absent),
                    ["figure_data/knn_nearest_target_%s.csv" % c for c in CONDS])
        return

    # Missing full_model_pu_scores.csv is a build-order violation, not an expected clean-clone state, and raises: pu_target_model.py produces it earlier in make tables.
    if not os.path.isfile(SCORES):
        raise SystemExit(f"{SCORES}\n  comes from `make tables` -- run that first.")

    os.makedirs(OUTDIR, exist_ok=True)
    reference, canon = anchor_set()
    var_name = np.loadtxt(os.path.join(KNN_DATA, "var_gene_names.txt"), dtype=str)

    for cond in CONDS:
        df = score_condition(cond, reference, canon, var_name)
        out = os.path.join(OUTDIR, f"knn_nearest_target_{cond}.csv")
        df.to_csv(out, index=False)
        s = df[df.gene == "STAT4"]
        extra = ""
        if len(s):
            r = s.iloc[0]
            extra = f" | STAT4 nearest={r.nearest_immune_target} cos={r.kNN_max_score:.3f}"
        print(f"[{cond}] {len(df)} genes, {int(df.is_approved_target.sum())} anchors{extra}")

    print(f"wrote {len(CONDS)} tables to {os.path.relpath(OUTDIR, ROOT)}/")
    print("next: add_knn_approved_drugs.py  (annotates these into figure_data/)")


if __name__ == "__main__":
    main()
