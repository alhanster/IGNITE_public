#!/usr/bin/env python3
"""Builds GPS per-gene inputs from Supplementary Table 16 of the GPS paper.

Writes:
  data/specificity/gps_per_gene.csv           gene, gene_id, gps_max_overall, gps_max_immune
  data/specificity/gps_drug_labeled_genes.csv the 453 genes GPS trained on
"""
import hashlib
import os
import sys

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

XLSX_PATH = os.path.join(ROOT, "data", "raw", "gps", "GPS_allgenes.xlsx")
PHECODE_PATH = os.path.join(ROOT, "data", "reference", "immune_phecodes_curated.csv")
UNIVERSE_PATH = os.path.join(ROOT, "data", "full_gene_list.with_perturbseq.tsv")

OUT_PER_GENE = os.path.join(ROOT, "data", "specificity", "gps_per_gene.csv")
OUT_LABELED = os.path.join(ROOT, "data", "specificity", "gps_drug_labeled_genes.csv")

SHEET = "Supplementary Table 16"
GENE_COL = "Gene"
SCORE_COL = "Genetic priority score (GPS)"
PHECODE_DESC_COL = "Phecode description"
PHECODE_INT_COL = "Phecode Integer"
INDICATION_COL = "Open Targets Indication"

# Table shape is asserted so a different release or partial download fails here, not downstream.
EXPECT_SRC_ROWS = 104870
EXPECT_SRC_GENES = 14899
EXPECT_SRC_PHECODES = 399

EXPECT_UNIVERSE = 19502
EXPECT_OVERALL_SCORED = 14677
EXPECT_LABELED = 453
EXPECT_IMMUNE_SCORED = 2539

# sha256 digests over the gene-sorted column at %.12f ("" for NaN); a failed digest means the source table or a derivation rule changed, recheck the 4.5a numbers if gps_max_overall moved.
SHA_OVERALL = "cc5f50ca84877a95cfbcefeb532f2e26e0207251fb755d659f7c0f7a16505db0"
SHA_IMMUNE = "c0f54f0bf965dbf28675ceeb18b3efd2cccce6dc2321ff927b03ce7e015fcd4e"
SHA_LABELED = "f164d8364ebb27d1e3df787933f101c912334fa209d429042771067ef11d921d"


def _digest_values(series):
    """sha256 of a float column computed at fixed precision, independent of float formatting."""
    return hashlib.sha256(
        ",".join("" if pd.isna(v) else "%.12f" % v for v in series).encode()
    ).hexdigest()


def _digest_genes(genes):
    return hashlib.sha256(",".join(sorted(genes)).encode()).hexdigest()


def main():
    if not os.path.isfile(XLSX_PATH):
        raise FileNotFoundError(
            "%s\n\n"
            "This is the published GPS supplementary table (Supplementary Table 16), the\n"
            "source for both gps_per_gene.csv and gps_drug_labeled_genes.csv. It is ~5.6 MB\n"
            "and gitignored under /data/raw/*, so a clean clone will not have it.\n"
            "Place the file at that path and re-run.\n" % XLSX_PATH
        )
    for p in (PHECODE_PATH, UNIVERSE_PATH):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

    src = pd.read_excel(XLSX_PATH, sheet_name=SHEET)
    for col in (GENE_COL, SCORE_COL, PHECODE_DESC_COL, PHECODE_INT_COL, INDICATION_COL):
        assert col in src.columns, "source table is missing column %r" % col
    assert len(src) == EXPECT_SRC_ROWS, \
        "source rows %d != expected %d -- wrong sheet or a re-released table" % (len(src), EXPECT_SRC_ROWS)
    assert src[GENE_COL].nunique() == EXPECT_SRC_GENES, \
        "source genes %d != expected %d" % (src[GENE_COL].nunique(), EXPECT_SRC_GENES)
    assert src[PHECODE_INT_COL].nunique() == EXPECT_SRC_PHECODES, \
        "source phecodes %d != expected %d" % (src[PHECODE_INT_COL].nunique(), EXPECT_SRC_PHECODES)
    assert src[SCORE_COL].notna().all(), "unexpected null GPS scores in the source table"

    immune_phecodes = pd.read_csv(PHECODE_PATH)["immune_phecode"].tolist()
    unknown = sorted(set(immune_phecodes) - set(src[PHECODE_DESC_COL].unique()))
    assert not unknown, \
        ("curated immune phecodes absent from the source table: %s\n"
         "A silent typo here would quietly shrink the immune filter." % unknown)

    universe = pd.read_csv(UNIVERSE_PATH, sep="\t", usecols=["gene", "gene_id"])
    assert universe.gene.is_unique, "gene universe must be unique to map onto"
    assert len(universe) == EXPECT_UNIVERSE, \
        "gene universe %d != expected %d" % (len(universe), EXPECT_UNIVERSE)

    max_overall = src.groupby(GENE_COL)[SCORE_COL].max()
    immune_src = src[src[PHECODE_DESC_COL].isin(immune_phecodes)]
    max_immune = immune_src.groupby(GENE_COL)[SCORE_COL].max()
    labeled = sorted(src.loc[src[INDICATION_COL].notna(), GENE_COL].unique())

    per_gene = universe.copy()
    per_gene["gps_max_overall"] = per_gene.gene.map(max_overall)
    per_gene["gps_max_immune"] = per_gene.gene.map(max_immune)

    n_overall = int(per_gene.gps_max_overall.notna().sum())
    n_immune = int(per_gene.gps_max_immune.notna().sum())
    assert n_overall == EXPECT_OVERALL_SCORED, \
        "genes with gps_max_overall %d != expected %d" % (n_overall, EXPECT_OVERALL_SCORED)
    assert n_immune == EXPECT_IMMUNE_SCORED, \
        "genes with gps_max_immune %d != expected %d" % (n_immune, EXPECT_IMMUNE_SCORED)
    assert len(labeled) == EXPECT_LABELED, \
        "drug-labeled genes %d != expected %d" % (len(labeled), EXPECT_LABELED)
    # gps_max_immune restricts gps_max_overall to a phecode subset, so it cannot exceed it.
    both = per_gene.dropna(subset=["gps_max_overall", "gps_max_immune"])
    assert (both.gps_max_immune <= both.gps_max_overall + 1e-9).all(), \
        "gps_max_immune exceeds gps_max_overall -- the immune filter is not a subset"

    ordered = per_gene.sort_values("gene").reset_index(drop=True)
    checks = [("gps_max_overall", _digest_values(ordered.gps_max_overall), SHA_OVERALL),
              ("gps_max_immune", _digest_values(ordered.gps_max_immune), SHA_IMMUNE),
              ("drug_labeled_genes", _digest_genes(labeled), SHA_LABELED)]
    bad = [(name, got, want) for name, got, want in checks if got != want]
    if bad:
        raise AssertionError(
            "derived values do not match the committed digests:\n"
            + "\n".join("  %-20s got %s\n  %-20s want %s" % (n, g, "", w) for n, g, w in bad)
            + "\nThe source table or a derivation rule changed. Establish which before\n"
              "updating a digest -- if gps_max_overall moved, every §4.5a number moves."
        )

    per_gene.to_csv(OUT_PER_GENE, index=False)
    pd.DataFrame({"gps_drug_labeled_gene": labeled}).to_csv(OUT_LABELED, index=False)

    print("  gps_per_gene.csv              %d genes | %d with gps_max_overall | %d with gps_max_immune"
          % (len(per_gene), n_overall, n_immune))
    print("  gps_drug_labeled_genes.csv    %d genes" % len(labeled))
    print("  all three value digests match the committed derivation")


if __name__ == "__main__":
    sys.exit(main())
