#!/usr/bin/env python
"""Stage-1 tables for the perturbational functional genomics figure (4 panels).

Builds the intermediate tables consumed by the figure renderer. See REPRODUCIBILITY.md,
Perturbational functional genomics gene universe, for the partition, denominator, and universe conventions.

Inputs  data/perturbseq/pu/pu_labels.csv
        data/full_gene_list.with_perturbseq.tsv
        data/perturbseq/sources/cytokine_significant_counts.csv
        data/perturbseq/sources/cytokine_receptor_significant_counts.csv
Outputs figure_data/th_de_long.csv                     gene, subset, neglog10_adjp, group
        figure_data/th_de_group_n.csv                  group, n
        figure_data/regulator_residual.csv             gene, group, residual
        figure_data/regulator_residual_group_n.csv     group, n
        figure_data/cytokine_counts.csv                gene, condition, n_sig, group
        figure_data/cytokine_group_n.csv               group, n
        figure_data/cytokine_receptor_counts.csv       gene, condition, n_sig, group
        figure_data/cytokine_receptor_group_n.csv      group, n
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
DATA = os.path.join(ROOT, "data")
FIGDATA = os.path.join(ROOT, "figure_data")

PU_LABELS = os.path.join(DATA, "perturbseq", "pu", "pu_labels.csv")
FEATURES  = os.path.join(DATA, "full_gene_list.with_perturbseq.tsv")
CYTOKINE  = os.path.join(DATA, "perturbseq", "sources",
                         "cytokine_significant_counts.csv")
RECEPTOR  = os.path.join(DATA, "perturbseq", "sources",
                         "cytokine_receptor_significant_counts.csv")

TARGET_LEVELS = ["approved", "in_trials", "other"]
TH_FEATURES = ["neglog10_adjp_Th1", "neglog10_adjp_Th2",
               "neglog10_adjp_Th17", "neglog10_adjp_Treg"]
TH_LABELS = ["Th1", "Th2", "Th17", "Treg"]
CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]


def drug_status_sets():
    """Builds status sets and the 19,502-gene universe from pu_role; kept separate to exclude out-of-universe genes."""
    lab = pd.read_csv(PU_LABELS)
    approved = set(lab.loc[lab["pu_role"] == "P", "gene"].unique())
    trial_only = set(lab.loc[lab["pu_role"] == "trial_heldout", "gene"].unique())
    universe = set(lab["gene"].unique())
    return approved, trial_only, universe


def assign_group(genes, approved, trial_only):
    return [
        "approved" if g in approved else "in_trials" if g in trial_only else "other"
        for g in genes
    ]


def counts_frame(series_of_groups):
    """Per-group counts in display order, as a tidy frame."""
    n = series_of_groups.value_counts()
    return pd.DataFrame([{"group": g, "n": int(n.get(g, 0))} for g in TARGET_LEVELS])


def restrict_to_universe(frame, universe, label, dropped):
    """Drops rows outside the universe and records the count of dropped rows into `dropped` for the log."""
    keep = frame["gene"].isin(universe)
    n_rows = int((~keep).sum())
    if n_rows:
        dropped.append((label, frame.loc[~keep, "gene"].nunique(), n_rows))
    return frame[keep].copy()


def build_cytokine_table(path, approved, trial_only, universe, label, dropped):
    """Shared builder for the panel c and panel d cytokine tables."""
    df = pd.read_csv(path)
    d = pd.DataFrame({
        "gene": df["gene"],
        "condition": df["culture_condition"],
        "n_sig": pd.to_numeric(df["n_sig"], errors="coerce"),
    })
    d = d[d["n_sig"].notna()].copy()
    d = restrict_to_universe(d, universe, label, dropped)
    # Kept as integer counts, not coerced via to_numeric (float64), so CSV output stays integral under pandas 3.0's trailing-.0 drop.
    d["n_sig"] = d["n_sig"].astype("Int64")
    d["group"] = assign_group(d["gene"], approved, trial_only)
    # One row per gene per condition, so the caption's n counts distinct genes, not rows.
    uniq = d.drop_duplicates(subset=["gene", "group"])
    return d[["gene", "condition", "n_sig", "group"]], counts_frame(uniq["group"])


def main():
    approved, trial_only, universe = drug_status_sets()
    fg = pd.read_csv(FEATURES, sep="\t")
    dropped = []   # (panel, n_genes, n_rows) for every out-of-universe drop

    blocks = []
    for col, label in zip(TH_FEATURES, TH_LABELS):
        blocks.append(pd.DataFrame({
            "gene": fg["gene"],
            "subset": label,
            "neglog10_adjp": pd.to_numeric(fg[col], errors="coerce"),
        }))
    th = pd.concat(blocks, ignore_index=True)
    th = th[th["neglog10_adjp"].notna()].copy()
    # A no-op today, since the feature TSV equals the universe, but applied to all four panels; if it ceases to be a no-op for only some, results diverge.
    th = restrict_to_universe(th, universe, "th_de", dropped)
    th["group"] = assign_group(th["gene"], approved, trial_only)
    # Counts distinct genes with a value in at least one subset, not the smaller non-NA n of any single violin.
    th_n = counts_frame(th.drop_duplicates(subset=["gene", "group"])["group"])

    rr = pd.DataFrame({
        "gene": fg["gene"],
        "residual": pd.to_numeric(fg["expected_n_regulators_residuals"],
                                  errors="coerce"),
    })
    rr = rr[rr["residual"].notna()].copy()
    rr = restrict_to_universe(rr, universe, "residual", dropped)
    rr["group"] = assign_group(rr["gene"], approved, trial_only)
    # Counts rows, not distinct genes; the two match here since there are no duplicate gene symbols.
    rr_n = counts_frame(rr["group"])

    cyt, cyt_n = build_cytokine_table(CYTOKINE, approved, trial_only, universe,
                                      "cytokine", dropped)
    rec, rec_n = build_cytokine_table(RECEPTOR, approved, trial_only, universe,
                                      "receptor", dropped)

    # Write outputs.
    os.makedirs(FIGDATA, exist_ok=True)
    th[["gene", "subset", "neglog10_adjp", "group"]].to_csv(
        os.path.join(FIGDATA, "th_de_long.csv"), index=False)
    th_n.to_csv(os.path.join(FIGDATA, "th_de_group_n.csv"), index=False)
    rr[["gene", "group", "residual"]].to_csv(
        os.path.join(FIGDATA, "regulator_residual.csv"), index=False)
    rr_n.to_csv(os.path.join(FIGDATA, "regulator_residual_group_n.csv"), index=False)
    cyt.to_csv(os.path.join(FIGDATA, "cytokine_counts.csv"), index=False)
    cyt_n.to_csv(os.path.join(FIGDATA, "cytokine_group_n.csv"), index=False)
    rec.to_csv(os.path.join(FIGDATA, "cytokine_receptor_counts.csv"), index=False)
    rec_n.to_csv(os.path.join(FIGDATA, "cytokine_receptor_group_n.csv"), index=False)

    print("perturbational functional-genomics tables written to figure_data/")
    for name, frame in [("th_de", th_n), ("residual", rr_n),
                        ("cytokine", cyt_n), ("receptor", rec_n)]:
        print("  %-10s %s" % (name, dict(zip(frame["group"], frame["n"]))))
    # The universe restriction must raise or log, never drop genes silently.
    if dropped:
        print("  out-of-universe rows dropped (absent from pu_labels.csv):")
        for label, n_genes, n_rows in dropped:
            print("    %-10s %d gene(s), %d row(s)" % (label, n_genes, n_rows))
    else:
        print("  out-of-universe rows dropped: none")


if __name__ == "__main__":
    main()
