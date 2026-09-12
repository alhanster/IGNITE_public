#!/usr/bin/env python3
"""build_vignette_panelB.py

Computes panel b of the STAT4 case vignette (fig_stat4_vignette.R) from the
committed feature matrix.

Panel b shows STAT4's percentile on four evidence axes, each computed against
measured (non-NaN) genes for that feature.

Input:  data/full_gene_list.with_perturbseq.tsv
Output: figure_data/vignette_panelB_axes.csv (radar panel, wrapped labels)
        figure_data/vignette_axes_supp.csv (Supplementary Table 14 rows)
Usage:  python src/analysis/build_vignette_panelB.py

See REPRODUCIBILITY.md, Vignette metadata, for the dropna
rationale and the two-file output rationale.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"

FEATURES = os.path.join(ROOT, "data", "full_gene_list.with_perturbseq.tsv")
OUT      = os.path.join(ROOT, "figure_data", "vignette_panelB_axes.csv")
OUT_SUPP = os.path.join(ROOT, "figure_data", "vignette_axes_supp.csv")
GENE     = "STAT4"

# axis label -> (feature column, evidence class); evidence keys fig_stat4_vignette.R's colour vector.
AXES = [
    ("Immune GWAS\n(genetic)",                    "gwas_score",                      "genetic"),
    ("Missense constraint\n(genetic)",            "mis.z_score",                     "genetic"),
    ("Th17 DE signal\n(observational)",           "neglog10_adjp_Th17",              "observational"),
    ("Trans-regulatory convergence\n(perturbational)",
                                                  "expected_n_regulators_residuals", "perturbational"),
]


def measured_percentile(series, value):
    """Percentile of `value` among measured (non-NaN) genes for this feature."""
    s = pd.to_numeric(series, errors="coerce")
    measured = s.dropna()                       # <-- load-bearing: measured genes only
    return (measured < value).mean() * 100.0


def main():
    fg = pd.read_csv(FEATURES, sep="\t")
    row = fg.loc[fg["gene"] == GENE]
    assert len(row) == 1, f"expected exactly one {GENE} row, found {len(row)}"

    rows = []
    for order, (label, col, evidence) in enumerate(AXES, start=1):
        assert col in fg.columns, f"missing feature column: {col}"
        val = pd.to_numeric(row[col], errors="coerce").iloc[0]
        assert pd.notna(val), f"{GENE} has no value for {col}"
        pct = measured_percentile(fg[col], val)
        rows.append({"axis": label, "feature": col, "percentile": round(pct, 1),
                     "evidence": evidence, "order": order})

    out = pd.DataFrame(rows, columns=["axis", "percentile", "evidence", "order"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {os.path.relpath(OUT, ROOT)}")

    # Supplementary Table 14 sheet: same four rows as the panel, label flattened, feature column added.
    supp = pd.DataFrame(rows, columns=["order", "axis", "feature", "percentile", "evidence"])
    supp["axis"] = supp["axis"].str.replace("\n", " ", regex=False)
    assert not supp["axis"].str.contains("\n").any(), "a label still carries a newline"
    assert list(supp["feature"]) == [col for _, col, _ in AXES], \
        "the feature column must be the AXES source columns, in order"
    supp.to_csv(OUT_SUPP, index=False)
    print(f"wrote {os.path.relpath(OUT_SUPP, ROOT)}")
    for r in rows:
        print(f"  {r['order']}. {r['axis'].replace(chr(10),' '):48s} {r['percentile']:5.1f}")


if __name__ == "__main__":
    main()
