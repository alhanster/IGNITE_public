"""See REPRODUCIBILITY.md, Other significance tests."""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from _stats import holm  # noqa: E402

from pathlib import Path
_REPO         = Path(__file__).resolve().parents[3]
IEI_CSV       = _REPO / "data" / "reference" / "IEI_gene_list.csv"
APPROVED_TXT  = _REPO / "data" / "data_drug" / "approved_target_genes.txt"
ALL_DRUGS_CSV = _REPO / "data" / "data_drug" / "all_drugs_approved_and_in_trial_by_gene_ot.csv"
DRUG_CSV      = _REPO / "data" / "reference" / "druggable_genome_gene_list.csv"
# Same 19,502-gene model universe panels b/c use (build_genetic_data_tables.py).
FULL_GENE_TSV = _REPO / "data" / "full_gene_list.tsv"
_OUTPUTS      = _REPO / "outputs"
_TABLES       = _OUTPUTS / "tables"
_FIGDATA      = _REPO / "figure_data"


def stars(p):
    """See REPRODUCIBILITY.md, Multiple-testing correction, for forest-plot stars."""
    if pd.isna(p):
        return ""
    return ("****" if p < 1e-4 else "***" if p < .001 else "**" if p < .01
            else "*" if p < .05 else "ns")


def load_lists():
    iei = set(pd.read_csv(IEI_CSV)['Gene'].dropna().astype(str).str.strip())
    with open(APPROVED_TXT) as f:
        approved = set(l.strip() for l in f if l.strip())
    dg = pd.read_csv(ALL_DRUGS_CSV)
    dg['gene_target'] = dg['gene_target'].dropna().astype(str).str.strip()
    # See REPRODUCIBILITY.md, Other significance tests.
    trial_staged = set(
        g for g in dg.loc[dg['furthest_stage'].isin(
            ['Phase 1', 'Phase 2', 'Phase 3']), 'gene_target'] if g)
    druggable = set(pd.read_csv(DRUG_CSV)['gene'].dropna().astype(str).str.strip())
    return iei, approved, trial_staged, druggable


def model_universe():
    """The 19,502-gene model universe (same file panels b/c partition by)."""
    df = pd.read_csv(FULL_GENE_TSV, sep="\t")
    return set(df['gene'].dropna().astype(str).str.strip())


def enrich_vs_reference(iei, group, reference, group_label, background):
    """OR of IEI-membership in `group` vs the `reference` background (2x2 Fisher)."""
    a = len(iei & group)                 # IEI in group
    b = len(group) - a                   # non-IEI in group
    c = len(iei & reference)             # IEI in reference background
    d = len(reference) - c               # non-IEI in reference background
    OR, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    if min(a, b, c, d) > 0:              # Woolf 95% CI on the odds ratio
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        lo, hi = np.exp(np.log(OR) - 1.96*se), np.exp(np.log(OR) + 1.96*se)
    else:
        lo, hi = np.nan, np.nan
    return dict(group=group_label, background=background, n=len(group), n_iei=a,
                p_iei=a/len(group), OR=OR, OR_lo=lo, OR_hi=hi, p=p, is_reference=False)


def build_enrichment_table():
    iei, approved, trial_staged, druggable = load_lists()

    trial_only = trial_staged - approved        # Phase 1/2/3 only; Unknown -> non-target
    all_pc     = model_universe()               # 19,502-gene model universe, no further union

    # Each background excludes both drug groups; only the reference set changes, and groups are never intersected with it.
    ref_all_genes = all_pc    - approved - trial_only
    ref_druggable = druggable - approved - trial_only

    groups = [("Targets of approved immune drugs", approved),
              ("Targets of immune drugs in trials", trial_only)]
    backgrounds = [("all_genes", ref_all_genes), ("druggable", ref_druggable)]

    rows = []
    for glabel, g in groups:
        for blabel, ref in backgrounds:
            rows.append(enrich_vs_reference(iei, g, ref, glabel, blabel))
    # One reference row per background, anchored at OR = 1 by definition, not a computed estimate.
    for glabel, ref, bg in [("All other genes", ref_all_genes, "all_genes"),
                            ("Druggable genome", ref_druggable, "druggable")]:
        rows.append(dict(group=glabel, background=bg,
                         n=len(ref), n_iei=len(iei & ref),
                         p_iei=len(iei & ref)/len(ref),
                         OR=1.0, OR_lo=np.nan, OR_hi=np.nan, p=np.nan, is_reference=True))
    res = pd.DataFrame(rows)

    # Readable companion
    bg_name = {"all_genes": "vs. all other genes",
               "druggable": "vs. druggable genome", "reference": "reference"}

    def or_txt(r):
        if r.is_reference:
            return "1.0 (reference)"
        return f"OR {r.OR:.1f} ({r.OR_lo:.1f}–{r.OR_hi:.1f})"
    readable = pd.DataFrame({
        'Group': res['group'],
        'Background': [bg_name[b] for b in res['background']],
        'Genes': res['n'],
        'IEI genes': [f"{r.n_iei} / {r.n} = {r.p_iei:.1%}" for r in res.itertuples()],
        'IEI enrichment (95% CI)': [or_txt(r) for r in res.itertuples()],
        'Fisher p': ["—" if r.is_reference else f"{r.p:.0e}" for r in res.itertuples()],
    })
    return res, readable


def write_markdown_table(res, path):
    """Write the wide two-background enrichment summary as a markdown table."""
    def or_cell(r):
        return f"{r.OR:.1f} ({r.OR_lo:.1f}–{r.OR_hi:.1f})"

    def row(r):  # a drug-group row: itertuples for both backgrounds
        all_r = r[r.background == "all_genes"].iloc[0]
        dg_r  = r[r.background == "druggable"].iloc[0]
        return (f"| {all_r.group} | {all_r.n:,} | {all_r.n_iei} | {all_r.p_iei:.1%} "
                f"| {or_cell(all_r)} | {or_cell(dg_r)} |")

    drug = res[~res.is_reference]
    lines = [
        "# IEI-gene enrichment among immune-drug targets",
        "",
        "| Group | n | IEI genes | IEI rate | OR vs all other genes (95% CI) | OR vs druggable genome (95% CI) |",
        "|---|---|---|---|---|---|",
    ]
    label = {"Targets of approved immune drugs": "Targets of approved immune drugs",
             "Targets of immune drugs in trials": "Targets of immune drugs in trials (not approved)"}
    for g in ["Targets of approved immune drugs", "Targets of immune drugs in trials"]:
        sub = drug[drug.group == g].copy()
        sub["group"] = label[g]
        lines.append(row(sub))
    # Reference rows (one per background), each anchored at OR = 1 in its own column by definition.
    ref_all = res[(res.is_reference) & (res.background == "all_genes")].iloc[0]
    ref_dg  = res[(res.is_reference) & (res.background == "druggable")].iloc[0]
    lines.append(f"| All other genes (reference) | {ref_all.n:,} | {ref_all.n_iei} "
                 f"| {ref_all.p_iei:.1%} | 1.0 (ref) | — |")
    lines.append(f"| Druggable genome (reference) | {ref_dg.n:,} | {ref_dg.n_iei} "
                 f"| {ref_dg.p_iei:.1%} | — | 1.0 (ref) |")
    lines += [
        "",
        "Enrichment of IEI genes among immune-drug targets. Each background excludes the two "
        "drug-target groups. Two-sided Fisher exact test; Woolf 95% CI.",
        "",
    ]
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    res, readable = build_enrichment_table()

    # ---- figure_data/: what the forest panel reads -------------------------------
    # Two columns the renderer consumes: `stars` (see stars() above) and `label_text`, the
    # finished row label. The renderer performs no number formatting itself.
    forest = res.copy()
    # Holm applies only to the four tested cells; the two OR=1 reference rows carry no p and are excluded rather than passed as NaN, which would silently take a real test's rank multiplier.
    tested = forest["p"].notna()
    assert int(tested.sum()) == 4, \
        f"the IEI family is four (group x background) cells, got {int(tested.sum())}"
    adj = pd.Series(np.nan, index=forest.index)
    adj[tested] = holm(forest.loc[tested, "p"].tolist())
    forest.insert(forest.columns.get_loc("p") + 1, "p_holm", adj)
    forest["stars"] = forest["p_holm"].map(stars)
    _REF_NAME = {"all_genes": "all other genes", "druggable": "druggable genome"}
    forest["label_text"] = [
        ("1.0 — %s (n=%s)" % (_REF_NAME[r.background], format(int(r.n), ",")))
        if bool(r.is_reference)
        else ("%.1f (%.1f–%.1f) %s" % (r.OR, r.OR_lo, r.OR_hi, r.stars))
        for r in forest.itertuples()
    ]
    _FIGDATA.mkdir(parents=True, exist_ok=True)
    forest.to_csv(_FIGDATA / "iei_enrichment_forest.csv", index=False)

    # Per-group display count uses max(n) across backgrounds, matching the renderer's tapply(res$n, res$group, max); a no-op since n matches across backgrounds.
    _CKEY = {"Targets of approved immune drugs": "approved",
             "Targets of immune drugs in trials": "in_trials"}
    grp_n = res.groupby("group")["n"].max()
    axis = pd.DataFrame(
        [{"ckey": ckey, "n": int(grp_n[name])} for name, ckey in _CKEY.items()],
        columns=["ckey", "n"])
    axis.to_csv(_FIGDATA / "iei_enrichment_axis.csv", index=False)

    # ---- outputs/: scratch. Neither of these reaches a figure. -------------------
    _TABLES.mkdir(parents=True, exist_ok=True)
    readable.to_csv(_TABLES / "iei_by_drug_status_readable.csv", index=False)
    write_markdown_table(res, _TABLES / "iei_enrichment_table.md")
    print(readable.to_string(index=False))
