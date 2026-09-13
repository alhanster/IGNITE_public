#!/usr/bin/env python3
"""Builds the per-panel univariate feature distribution tables behind Figures 2 and 3, plus IEI enrichment and per-panel group sizes.

Runs as part of `make tables`, after build_ranked_atlas_table.py; reshapes outputs already written by make_enrichment_table.py, build_genetic_data_tables.py and build_perturbational_fg_tables.py. Reshaping only, with six exceptions: univariate_gwas_tests.csv, univariate_mis_z_tests.csv, univariate_th_de_fdr_enrichment.csv, univariate_residual_tests.csv, univariate_cytokine_distribution.csv and univariate_cytokine_regulator_fraction.csv test hypotheses no upstream script asks. All six are computed here rather than upstream because their inputs are the panels' own committed sheets, so none can disagree with Figures 2 and 3 as drawn, and all rebuild on a clone with no perturb-seq inputs. Between them they carry the statistics for every panel of both figures except 2a, whose Fisher enrichment make_enrichment_table.py already writes.

Outputs (15 sheets): univariate_iei_enrichment.csv, univariate_gwas.csv, univariate_gwas_pct_gt0.csv, univariate_gwas_tests.csv, univariate_mis_z.csv, univariate_mis_z_tests.csv, univariate_th_de.csv, univariate_th_de_fdr_enrichment.csv, univariate_residual.csv, univariate_residual_tests.csv, univariate_cytokines.csv, univariate_cytokine_receptors.csv, univariate_cytokine_distribution.csv, univariate_cytokine_regulator_fraction.csv, univariate_group_n.csv.
"""
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, kruskal, mannwhitneyu, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from _stats import holm  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"
FIGDATA = os.path.join(ROOT, "figure_data")

# Maps the Fig 2/3 group-label spelling to the spelling these output sheets publish.
REMAP = {"approved": "approved", "in_trials": "in-trial", "other": "non-target"}
# ranked_atlas.csv already uses this spelling; kept as an identity map, the seam regroup() asserts through.
REMAP_ATLAS = {"approved": "approved", "in-trial": "in-trial", "non-target": "non-target"}
# Six rows in iei_enrichment_forest.csv; two are BACKGROUND reference rows, not drug-status groups, hence four entries.
REMAP_ENRICH = {"Targets of approved immune drugs": "approved immune targets",
                "Targets of immune drugs in trials": "in-trial immune targets",
                "Druggable genome": "druggable genome",
                "All other genes": "all other genes"}


# Fig 3a's own threshold, applied to the panel's y-value: neglog10_adjp >= -log10(0.10) == 1.0.
# No row of th_de_long.csv sits exactly on the boundary today, so >= and > agree; >= is the
# inclusive reading of "reaching 10% FDR" and is the one the caption claims.
FDR_ALPHA = 0.10
FDR_CUT = -math.log10(FDR_ALPHA)
TH_LABELS = ["Th1", "Th2", "Th17", "Treg"]
# Each target class is its own pre-specified family of four subsets, corrected independently:
# the claim is per class ("approved targets are enriched in Th1"), so a class's four subsets are
# what must hold together, and the two classes are not a single screen over eight cells. See
# REPRODUCIBILITY.md.
FDR_TARGET_CLASSES = ["approved", "in-trial"]


def fd(name):
    return os.path.join(FIGDATA, name)


def woolf_ci(OR, a, b, c, d):
    """Woolf 95% CI on a 2x2 odds ratio, spelled as make_enrichment_table.py spells it.

    Both Fisher sheets in this module call this rather than inlining the formula: the algebraically
    equal `OR * exp(+-z*se)` disagrees with `exp(log(OR) +- z*se)` in the last bit on roughly half
    of these cells, so two spellings of one interval would leave figure_data/ unable to survive a
    later tidy-up that merged them -- `make verify-tables` would fail on a last digit with nothing
    naming the cause.

    The 1.96 is the repo-wide literal, matching make_enrichment_table.py, build_temporal_holdout.py
    and compute_delong_specificity.py; it is deliberately NOT derived from CI_LEVEL, which governs
    only the rank-based intervals below. A zero cell leaves the interval undefined, as upstream.
    """
    if min(a, b, c, d) <= 0:
        return float("nan"), float("nan")
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return math.exp(math.log(OR) - 1.96 * se), math.exp(math.log(OR) + 1.96 * se)


# Fig 3c/3d: the stimulation series, in the order the panels facet it, and the two target classes
# each tested against non-target. Each condition is one pre-specified question -- "at this
# timepoint, is a target class more likely to regulate anything?" -- so each condition's two
# tests are one Holm family. That boundary was fixed by the panel's facet structure before the
# tests were run, not chosen after seeing which cell survived.
CYT_CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]
CYT_TARGET_CLASSES = ["approved", "in-trial"]
CYT_REFERENCE = "non-target"
CYT_PANELS = [("3c cytokines", "univariate_cytokines.csv"),
              ("3d cytokine receptors", "univariate_cytokine_receptors.csv")]


def check_cyt_labels(panel, df):
    """Names a renamed condition or group, as the Fig 3a and 3b builders do for theirs.

    Without it a label that stopped matching REMAP surfaces as an empty-Series max() or a
    zero denominator several frames down, with nothing naming which label moved.
    """
    for col, want in (("condition", CYT_CONDITIONS),
                      ("group", CYT_TARGET_CLASSES + [CYT_REFERENCE])):
        missing = set(want) - set(df[col])
        assert not missing, f"{panel}: sheet is missing {col} label(s): {sorted(missing)}"


def cytokine_distribution(sheets):
    """Per condition x group shape of the Fig 3c/3d counts.

    Exists because the Fig 3b rank-test template is degenerate on these data: at 88-94% zeros
    every median, Hodges-Lehmann shift and interval bound is exactly 0, so a location-test sheet
    would carry four columns of zeros and say nothing. What does separate the groups is how often
    a gene regulates anything at all, and how heavy the tail is when it does -- this sheet reports
    both, and `n_ge3` and `max` are what let a reader see that panel C/D's mean diamond rides on a
    handful of genes.
    """
    rows = []
    for panel, fname in CYT_PANELS:
        df = sheets[fname]
        check_cyt_labels(panel, df)
        for cond in CYT_CONDITIONS:
            for grp in CYT_TARGET_CLASSES + [CYT_REFERENCE]:
                s = df.loc[(df["condition"] == cond) & (df["group"] == grp), "n_sig"]
                nz = s[s > 0]
                rows.append({"panel": panel, "condition": cond, "group": grp,
                             "n": len(s), "n_reg": int((s > 0).sum()),
                             "pct_reg": 100 * float((s > 0).mean()),
                             "pct_zero": 100 * float((s == 0).mean()),
                             "mean": float(s.mean()),
                             "mean_nonzero": float(nz.mean()) if len(nz) else "",
                             "max": int(s.max()), "n_ge3": int((s >= 3).sum())})
    return pd.DataFrame(rows)


def cytokine_regulator_fraction(sheets):
    """Fisher tests on the fraction of genes regulating at least one target, per condition.

    Dichotomizes `n_sig` at zero and asks, within a condition, whether a target class regulates
    anything more often than non-target genes. This is the first half of a two-part reading of
    panels C and D and the part these distributions can actually support; the magnitude half is
    left to `mean_nonzero` on the distribution sheet rather than tested, since the non-zero
    subsets run to a few dozen genes.

    Holm spans the two target classes within one condition -- six families of two.
    """
    rows = []
    for panel, fname in CYT_PANELS:
        df = sheets[fname]
        check_cyt_labels(panel, df)
        for cond in CYT_CONDITIONS:
            dc = df[df["condition"] == cond]
            ref = dc.loc[dc["group"] == CYT_REFERENCE, "n_sig"]
            c, d = int((ref > 0).sum()), int((ref == 0).sum())
            family = []
            for grp in CYT_TARGET_CLASSES:
                s = dc.loc[dc["group"] == grp, "n_sig"]
                a, b = int((s > 0).sum()), int((s == 0).sum())
                OR, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
                lo, hi = woolf_ci(float(OR), a, b, c, d)
                family.append({"panel": panel, "condition": cond, "group": grp,
                               "n": a + b, "n_reg": a, "pct_reg": 100 * a / (a + b),
                               "n_ref": c + d, "n_reg_ref": c, "pct_reg_ref": 100 * c / (c + d),
                               "OR": float(OR), "OR_lo": lo, "OR_hi": hi,
                               "p": float(p)})
            assert len(family) == len(CYT_TARGET_CLASSES), \
                f"{panel}/{cond}: Holm family has {len(family)} members, expected {len(CYT_TARGET_CLASSES)}"
            for row, p_holm in zip(family, holm([r["p"] for r in family])):
                row["p_holm"] = p_holm
                row["stars"] = stars(p_holm)
            rows += family
    return pd.DataFrame(rows)


# Fig 3b pairwise comparisons, in the order the panel draws its violins. The reference group is
# second in each pair, so a positive Hodges-Lehmann shift means the first group carries the
# larger residual burden.
RESIDUAL_PAIRS = [("approved", "non-target"),
                  ("in-trial", "non-target"),
                  ("approved", "in-trial")]
CI_LEVEL = 0.95


def stars(p):
    """The repo's significance ladder, as build_score_by_group.py and build_leakage_controlled.py spell it."""
    return ("****" if p < 1e-4 else "***" if p < 1e-3
            else "**" if p < 1e-2 else "*" if p < 0.05 else "ns")


def hodges_lehmann(x, y, conf=CI_LEVEL):
    """HL shift estimate and its Wilcoxon rank-sum (Moses) confidence interval.

    The estimate is the median of all n1*n2 pairwise differences; the interval is the pair of
    order statistics of those differences cut off by the rank-sum null's critical value. Chosen
    over a bootstrap deliberately: figure_data/ is verified byte-identical on every rebuild, and
    an interval that depends on a seed and a draw count cannot satisfy that. It agrees with R's
    wilcox.test(conf.int = TRUE, exact = FALSE) to five decimal places but is not identical to it:
    that call root-finds a continuous bound, where this selects an order statistic.
    """
    n1, n2 = len(x), len(y)
    d = np.sort((np.asarray(x)[:, None] - np.asarray(y)[None, :]).ravel())
    k = d.size
    z = norm.ppf(0.5 + conf / 2)
    c = int(np.floor(n1 * n2 / 2 - z * math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    c = min(max(c, 0), k // 2 - 1)      # clamp: a group too small to bound returns the extremes
    return float(np.median(d)), float(d[c]), float(d[k - 1 - c])


def prob_superiority(x, y, conf=CI_LEVEL):
    """P(sup) = U / (n1*n2) with its Hanley-McNeil normal-approximation interval.

    The rank-biserial correlation is deliberately not emitted beside this: it equals
    2*P(sup) - 1 identically, so the two columns would carry one number twice.
    """
    n1, n2 = len(x), len(y)
    U = float(mannwhitneyu(x, y, alternative="two-sided").statistic)
    A = U / (n1 * n2)
    q1 = A / (2 - A)
    q2 = 2 * A ** 2 / (1 + A)
    se = math.sqrt((A * (1 - A) + (n1 - 1) * (q1 - A ** 2) + (n2 - 1) * (q2 - A ** 2)) / (n1 * n2))
    z = norm.ppf(0.5 + conf / 2)
    return U, A, max(0.0, A - z * se), min(1.0, A + z * se)


def residual_group_tests(res):
    """Kruskal-Wallis omnibus plus three pairwise Mann-Whitney tests on the Fig 3b residual.

    `res` is the published-vocabulary residual sheet (approved / in-trial / non-target). Long
    format, omnibus first, matching score_by_group_tests.csv: the omnibus lives in the sheet as
    a row rather than in a caption, so the sheet stands alone once it is an xlsx tab. The
    omnibus leaves every pairwise column empty and is left uncorrected, as that sheet's
    Kruskal-Wallis and ANOVA rows are.

    Holm spans the three pairwise tests only -- one family, matching build_score_by_group.py's
    identically shaped Figure 7 comparison. See REPRODUCIBILITY.md.
    """
    g = {k: d["residual"].to_numpy() for k, d in res.groupby("group")}
    missing = {a for pair in RESIDUAL_PAIRS for a in pair} - set(g)
    assert not missing, f"residual sheet is missing group(s): {sorted(missing)}"

    H, p_kw = kruskal(*(g[k] for k in ("approved", "in-trial", "non-target")))
    rows = [{"test": "Kruskal-Wallis", "group_a": "", "group_b": "", "n_a": "", "n_b": "",
             "median_a": "", "median_b": "", "hl_shift": "", "hl_lo": "", "hl_hi": "",
             "p_sup": "", "p_sup_lo": "", "p_sup_hi": "", "statistic": "H",
             "value": float(H), "p": float(p_kw), "p_holm": "", "stars": ""}]

    pairwise = []
    for a, b in RESIDUAL_PAIRS:
        x, y = g[a], g[b]
        hl, hl_lo, hl_hi = hodges_lehmann(x, y)
        U, A, a_lo, a_hi = prob_superiority(x, y)
        p = float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
        pairwise.append({"test": "Mann-Whitney (two-sided)", "group_a": a, "group_b": b,
                         "n_a": len(x), "n_b": len(y),
                         "median_a": float(np.median(x)), "median_b": float(np.median(y)),
                         "hl_shift": hl, "hl_lo": hl_lo, "hl_hi": hl_hi,
                         "p_sup": A, "p_sup_lo": a_lo, "p_sup_hi": a_hi,
                         "statistic": "U", "value": U, "p": p})

    assert len(pairwise) == 3, f"Holm family has {len(pairwise)} members, expected 3"
    for row, p_holm in zip(pairwise, holm([r["p"] for r in pairwise])):
        row["p_holm"] = p_holm
        row["stars"] = stars(p_holm)

    return pd.DataFrame(rows + pairwise)


# Fig 2b/2c: the same three drug-status comparisons as Fig 3b, on the two genetic-prior panels.
GWAS_VALUE = "gwas_score"
MIS_Z_VALUE = "mis_z_score"
GROUP_ORDER = ["approved", "in-trial", "non-target"]

# Declared rather than inherited from whatever order the row dicts happen to be built in, as
# th_de_fdr_enrichment declares its own: these sheets are read beside univariate_residual_tests.csv
# and carry its column order with `genes` second, Fig 2b's presence-only columns sitting after n_b.
# PROVENANCE.md quotes these headers verbatim, so a reordering here has to be made there too.
RANK_COLUMNS = ["test", "genes", "group_a", "group_b", "n_a", "n_b",
                "median_a", "median_b", "hl_shift", "hl_lo", "hl_hi",
                "p_sup", "p_sup_lo", "p_sup_hi", "statistic", "value", "p", "p_holm", "stars"]
GWAS_COLUMNS = (RANK_COLUMNS[:6]
                + ["n_gt0_a", "n_gt0_b", "pct_gt0_a", "pct_gt0_b", "OR", "OR_lo", "OR_hi"]
                + RANK_COLUMNS[6:])


def frame(rows, columns):
    """Rows to a DataFrame, filling every column a row does not carry with an empty string.

    Empty rather than absent, for the reason residual_group_tests spells its omnibus row out in
    full: a missing key becomes NaN, which makes the whole column float, and `n_a` would then be
    written `340.0` where that sibling sheet writes `206`. package_supplementary_tables.py copies
    these files with dtype=str, so the CSV's formatting is what reaches the workbook.
    """
    return pd.DataFrame([{c: r.get(c, "") for c in columns} for r in rows], columns=columns)


def kw_pairwise_block(g, genes_label):
    """Kruskal-Wallis omnibus row plus Holm-corrected pairwise Mann-Whitney rows.

    `g` maps group label to a value array. Row shape follows univariate_residual_tests.csv, so
    every rank-test sheet in this table reads the same way; `genes` records which genes the block
    was computed on, because Figure 2b's two blocks are computed on different subsets. The omnibus
    is a row rather than a caption and is left uncorrected, as in score_by_group_tests.csv.
    """
    missing = set(GROUP_ORDER) - set(g)
    assert not missing, f"{genes_label}: sheet is missing group(s): {sorted(missing)}"

    H, p_kw = kruskal(*(g[k] for k in GROUP_ORDER))
    rows = [{"test": "Kruskal-Wallis", "genes": genes_label, "group_a": "", "group_b": "",
             "statistic": "H", "value": float(H), "p": float(p_kw)}]
    family = []
    for a, b in RESIDUAL_PAIRS:
        x, y = g[a], g[b]
        hl, hl_lo, hl_hi = hodges_lehmann(x, y)
        U, A, a_lo, a_hi = prob_superiority(x, y)
        family.append({"test": "Mann-Whitney (two-sided)", "genes": genes_label,
                       "group_a": a, "group_b": b, "n_a": len(x), "n_b": len(y),
                       "median_a": float(np.median(x)), "median_b": float(np.median(y)),
                       "hl_shift": hl, "hl_lo": hl_lo, "hl_hi": hl_hi,
                       "p_sup": A, "p_sup_lo": a_lo, "p_sup_hi": a_hi,
                       "statistic": "U", "value": U,
                       "p": float(mannwhitneyu(x, y, alternative="two-sided").pvalue)})
    assert len(family) == len(RESIDUAL_PAIRS), \
        f"{genes_label}: Holm family has {len(family)} members, expected {len(RESIDUAL_PAIRS)}"
    for row, p_holm in zip(family, holm([r["p"] for r in family])):
        row["p_holm"] = p_holm
        row["stars"] = stars(p_holm)
    return rows + family


def gwas_tests(sheets):
    """Fig 2b in two parts: who carries any immune-GWAS signal, and how much when they do.

    The panel's y-axis reads "genes with signal (>0)", so a rank test describing the panel as
    drawn has to be restricted the same way. But `gwas_score` is 100% non-missing and not imputed
    (see figure_data/feature_dictionary.csv), so a zero is a measured absence of association, not
    a gap -- the dropped genes are data, and the fraction above zero is itself testable. Hence two
    blocks: Fisher on all genes for the presence question, rank tests on the subset for magnitude.

    The two blocks answer different questions and are corrected separately, three tests each.
    Reading only one of them is what misleads: on presence the two target classes are
    indistinguishable from each other, while on magnitude approved exceeds in-trial. See
    REPRODUCIBILITY.md.
    """
    df = sheets["univariate_gwas.csv"]
    g_all = {k: d[GWAS_VALUE].to_numpy() for k, d in df.groupby("group")}
    missing = set(GROUP_ORDER) - set(g_all)
    assert not missing, f"gwas sheet is missing group(s): {sorted(missing)}"

    presence = []
    for a, b in RESIDUAL_PAIRS:
        x, y = g_all[a], g_all[b]
        ha, na = int((x > 0).sum()), len(x)
        hb, nb = int((y > 0).sum()), len(y)
        OR, p = fisher_exact([[ha, na - ha], [hb, nb - hb]], alternative="two-sided")
        OR_lo, OR_hi = woolf_ci(OR, ha, na - ha, hb, nb - hb)
        presence.append({"test": "Fisher exact (any signal >0)", "genes": "all genes",
                         "group_a": a, "group_b": b, "n_a": na, "n_b": nb,
                         "n_gt0_a": ha, "n_gt0_b": hb,
                         "pct_gt0_a": 100 * ha / na, "pct_gt0_b": 100 * hb / nb,
                         "OR": float(OR), "OR_lo": OR_lo, "OR_hi": OR_hi, "p": float(p)})
    assert len(presence) == len(RESIDUAL_PAIRS), \
        f"presence Holm family has {len(presence)} members, expected {len(RESIDUAL_PAIRS)}"
    for row, p_holm in zip(presence, holm([r["p"] for r in presence])):
        row["p_holm"] = p_holm
        row["stars"] = stars(p_holm)

    pos = df[df[GWAS_VALUE] > 0]
    g_pos = {k: d[GWAS_VALUE].to_numpy() for k, d in pos.groupby("group")}
    return frame(presence + kw_pairwise_block(g_pos, "with signal (>0)"), GWAS_COLUMNS)


def mis_z_tests(sheets):
    """Fig 2c: the same rank battery on the missense Z-score, with no subsetting.

    2b splits because `gwas_score` carries a 54-64% mass at exactly zero, not because of anything
    to do with missingness; `mis_z_score` has no zero mass at all (no gene scores exactly 0), so
    the whole panel is one block and every effect-size column is well defined. The feature is
    *not* the fully observed one -- it is 99.62% non-missing and median-imputed during fitting
    (figure_data/feature_dictionary.csv), the opposite of `gwas_score` on both counts. The 75
    genes without a score, all non-target, are dropped upstream by build_genetic_data_tables.py,
    which is why this sheet's non-target n is 18,360 against 2b's 18,435.
    """
    df = sheets["univariate_mis_z.csv"]
    g = {k: d[MIS_Z_VALUE].to_numpy() for k, d in df.groupby("group")}
    return frame(kw_pairwise_block(g, "all genes"), RANK_COLUMNS)


def th_de_fdr_enrichment(th):
    """Per-subset 2x2 Fisher tests of 10% FDR attainment, each target class vs non-targets.

    `th` is the published-vocabulary T-helper sheet (approved / in-trial / non-target). Four
    rows per subset, in display order: the marginal rate over all tested genes, the two target
    classes, then the non-target reference. The marginal and reference rows carry no odds ratio
    against themselves and leave the inferential columns empty, as univariate_iei_enrichment.csv
    does for its background rows.

    n_tested is per subset, not per panel: a gene absent from one subset's DE table is dropped
    from that subset only, so these counts are smaller than the facet-strip n and differ between
    subsets. Holm runs within each target class over its four subsets -- two families of four,
    not one of eight -- so p_holm for an approved row is never a function of an in-trial p-value.
    """
    hit = th["neglog10_adjp"] >= FDR_CUT
    rows, tests = [], []
    for subset in TH_LABELS:
        m = th["subset"] == subset
        assert m.any(), f"no rows for subset {subset}"
        ref = m & (th["group"] == "non-target")
        c, d = int(hit[ref].sum()), int((~hit[ref]).sum())

        rows.append({"subset": subset, "group": "all tested genes",
                     "n_tested": int(m.sum()), "n_hit": int(hit[m].sum())})
        for label in FDR_TARGET_CLASSES:
            g = m & (th["group"] == label)
            a, b = int(hit[g].sum()), int((~hit[g]).sum())
            OR, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            lo, hi = woolf_ci(OR, a, b, c, d)
            row = {"subset": subset, "group": label, "n_tested": a + b, "n_hit": a,
                   "OR": OR, "OR_lo": lo, "OR_hi": hi, "p": p}
            rows.append(row)
            tests.append(row)
        rows.append({"subset": subset, "group": "non-target",
                     "n_tested": c + d, "n_hit": c, "OR": 1.0})

    # One Holm family per target class. Asserted rather than assumed: a subset that stopped
    # contributing a row would otherwise silently shrink a family and inflate its adjusted p.
    for label in FDR_TARGET_CLASSES:
        family = [r for r in tests if r["group"] == label]
        assert len(family) == len(TH_LABELS), \
            f"{label}: Holm family has {len(family)} members, expected {len(TH_LABELS)}"
        for row, p_holm in zip(family, holm([r["p"] for r in family])):
            row["p_holm"] = p_holm

    out = pd.DataFrame(rows, columns=["subset", "group", "n_tested", "n_hit", "pct_hit",
                                      "OR", "OR_lo", "OR_hi", "p", "p_holm"])
    out["pct_hit"] = 100 * out["n_hit"] / out["n_tested"]
    return out


def regroup(df, mapping, col="group"):
    unknown = set(df[col]) - set(mapping)
    assert not unknown, f"unmapped group label(s): {sorted(unknown)}"
    out = df.copy()
    out[col] = out[col].map(mapping)
    return out


def main():
    sheets = {}

    # Fig 2a: IEI enrichment vs both backgrounds.
    iei = regroup(pd.read_csv(fd("iei_enrichment_forest.csv")), REMAP_ENRICH)
    sheets["univariate_iei_enrichment.csv"] = iei[
        ["group", "background", "n", "n_iei", "OR", "OR_lo", "OR_hi", "p", "p_holm"]]
    assert len(iei) == 6, f"expected 6 enrichment rows, got {len(iei)}"

    # Fig 2b: immune-GWAS score, unfiltered.
    atlas = pd.read_csv(fd("ranked_atlas.csv"), usecols=["gene", "drug_status", "gwas_score"])
    assert atlas["gwas_score"].notna().all(), "gwas_score is expected to be complete in the atlas"
    gwas = regroup(atlas.rename(columns={"drug_status": "group"}), REMAP_ATLAS)
    sheets["univariate_gwas.csv"] = gwas[["gene", "group", "gwas_score"]]
    assert len(gwas) == 19502, f"expected 19,502 genes, got {len(gwas)}"

    # pct_gt0 is stored as a fraction upstream; scaled by 100 here to match its name, not the zero-inflated figure quoted elsewhere.
    pct = regroup(pd.read_csv(fd("gwas_zero_inflation.csv")), REMAP)
    pct["pct_gt0"] = pct["pct_gt0"] * 100
    sheets["univariate_gwas_pct_gt0.csv"] = pct

    passthrough = {
        "univariate_mis_z.csv": ("mis_z_by_drug_status.csv", ["gene", "group", "mis_z_score"], 19427),
        "univariate_th_de.csv": ("th_de_long.csv", ["gene", "subset", "neglog10_adjp", "group"], 38827),
        "univariate_residual.csv": ("regulator_residual.csv", ["gene", "group", "residual"], 10268),
        "univariate_cytokines.csv": ("cytokine_counts.csv", ["gene", "condition", "n_sig", "group"], 33861),
        "univariate_cytokine_receptors.csv": ("cytokine_receptor_counts.csv",
                                              ["gene", "condition", "n_sig", "group"], 33861),
    }
    for out_name, (src, cols, n_expected) in passthrough.items():
        df = regroup(pd.read_csv(fd(src)), REMAP)
        assert len(df) == n_expected, f"{src}: expected {n_expected} rows, got {len(df)}"
        sheets[out_name] = df[cols]

    # Fig 3a: computed off the published sheet, so the tests describe the panel as drawn.
    fdr = th_de_fdr_enrichment(sheets["univariate_th_de.csv"])
    assert len(fdr) == 16, f"expected 16 rows (4 subsets x 4 groups), got {len(fdr)}"
    # The three group rows must partition each subset's tested genes; a group label that stopped
    # matching REMAP would otherwise leave a subset silently short of its own marginal.
    for subset in TH_LABELS:
        s = fdr[fdr["subset"] == subset].set_index("group")["n_tested"]
        parts = int(s[["approved", "in-trial", "non-target"]].sum())
        assert parts == int(s["all tested genes"]), \
            f"{subset}: groups sum to {parts}, marginal is {int(s['all tested genes'])}"
    sheets["univariate_th_de_fdr_enrichment.csv"] = fdr

    # Fig 2b/2c: tested off the published sheets, like the Figure 3 test sheets below.
    gtests = gwas_tests(sheets)
    ztests = mis_z_tests(sheets)
    assert len(gtests) == 7, f"expected 3 Fisher + 1 omnibus + 3 pairwise rows, got {len(gtests)}"
    assert len(ztests) == 4, f"expected 1 omnibus + 3 pairwise rows, got {len(ztests)}"
    # The presence rows must reproduce the published percentages exactly; if the two ever drifted,
    # Table 4 would carry a percentage and a test of that percentage that disagreed. Read from the
    # sheet this run writes, not from the file the last run left on disk: that file is rewritten
    # below, so reading it would check the new tests against a stale copy of their own sibling.
    pub = sheets["univariate_gwas_pct_gt0.csv"].set_index("group")["pct_gt0"]
    for _, r in gtests[gtests["test"].str.startswith("Fisher")].iterrows():
        for side in ("a", "b"):
            assert abs(r[f"pct_gt0_{side}"] - pub[r[f"group_{side}"]]) < 1e-9, \
                f"{r[f'group_{side}']}: pct above zero disagrees with univariate_gwas_pct_gt0.csv"
    # A median or shift of exactly 0 across a whole block is the signature of testing a
    # zero-inflated column unrestricted -- the failure mode that keeps 2b's rank tests on the
    # subset the panel draws. Guard it rather than rediscover it.
    for name, t in (("2b", gtests), ("2c", ztests)):
        mw = t[t["test"].str.startswith("Mann-Whitney")]
        assert not (mw["hl_shift"] == 0).all(), f"{name}: every Hodges-Lehmann shift is zero"
        assert not (mw["median_a"] == 0).all(), f"{name}: every median is zero"
    sheets["univariate_gwas_tests.csv"] = gtests
    sheets["univariate_mis_z_tests.csv"] = ztests

    # Fig 3b: same principle as above -- tested off the published sheet, not off data/.
    rtests = residual_group_tests(sheets["univariate_residual.csv"])
    assert len(rtests) == 4, f"expected 1 omnibus + 3 pairwise rows, got {len(rtests)}"
    # Each pairwise row's n and medians must agree with the source sheet's own groups, which
    # catches a pair written against the wrong group before the numbers reach a caption.
    src_res = sheets["univariate_residual.csv"]
    for _, r in rtests[rtests["test"] != "Kruskal-Wallis"].iterrows():
        for side, grp_name in (("a", r["group_a"]), ("b", r["group_b"])):
            col = src_res.loc[src_res["group"] == grp_name, "residual"]
            assert r[f"n_{side}"] == len(col), \
                f"{grp_name}: sheet has {len(col)} rows, test row says {r[f'n_{side}']}"
            assert abs(r[f"median_{side}"] - col.median()) < 1e-12, \
                f"{grp_name}: median disagrees with the source sheet"
    sheets["univariate_residual_tests.csv"] = rtests

    # Fig 3c/3d: same principle again -- both read the published panel sheets, not data/.
    cdist = cytokine_distribution(sheets)
    cfrac = cytokine_regulator_fraction(sheets)
    n_cells = len(CYT_PANELS) * len(CYT_CONDITIONS)
    assert len(cdist) == n_cells * (len(CYT_TARGET_CLASSES) + 1), \
        f"expected {n_cells * 3} distribution rows, got {len(cdist)}"
    assert len(cfrac) == n_cells * len(CYT_TARGET_CLASSES), \
        f"expected {n_cells * 2} Fisher rows, got {len(cfrac)}"
    # Each Fisher row's counts must match the distribution sheet's own cell for the same
    # panel/condition/group, so the two sheets can never drift into disagreeing with each other.
    key = ["panel", "condition", "group"]
    merged = cfrac.merge(cdist[key + ["n", "n_reg"]], on=key, suffixes=("", "_dist"))
    assert len(merged) == len(cfrac), "a Fisher row has no matching distribution cell"
    assert (merged["n"] == merged["n_dist"]).all() and (merged["n_reg"] == merged["n_reg_dist"]).all(), \
        "Fisher counts disagree with the distribution sheet"
    sheets["univariate_cytokine_distribution.csv"] = cdist
    sheets["univariate_cytokine_regulator_fraction.csv"] = cfrac

    # genetic_data_group_n.csv already carries a panel column; the four Fig 3 files are group,n with panel encoded only in the filename, so the label is attached here.
    parts = [regroup(pd.read_csv(fd("genetic_data_group_n.csv")), REMAP)[["panel", "group", "n"]]]
    for panel, src in [("th_de", "th_de_group_n.csv"),
                       ("residual", "regulator_residual_group_n.csv"),
                       ("cytokines", "cytokine_group_n.csv"),
                       ("cytokine_receptors", "cytokine_receptor_group_n.csv")]:
        part = regroup(pd.read_csv(fd(src)), REMAP)
        parts.append(part.assign(panel=panel)[["panel", "group", "n"]])
    group_n = pd.concat(parts, ignore_index=True)
    sheets["univariate_group_n.csv"] = group_n

    # Reconciles each panel's counts against its source sheet, which catches a mislabelled panel a row-count check alone would miss.
    for panel, sheet, distinct in [("mis_z", "univariate_mis_z.csv", True),
                                   ("th_de", "univariate_th_de.csv", True),
                                   ("residual", "univariate_residual.csv", False),
                                   ("cytokines", "univariate_cytokines.csv", True),
                                   ("cytokine_receptors", "univariate_cytokine_receptors.csv", True)]:
        stated = int(group_n.loc[group_n["panel"] == panel, "n"].sum())
        df = sheets[sheet]
        actual = df["gene"].nunique() if distinct else len(df)
        assert stated == actual, f"{panel}: group_n sums to {stated}, sheet has {actual}"

    for name, df in sheets.items():
        df.to_csv(fd(name), index=False)
        print(f"wrote figure_data/{name}  ({len(df):,} rows x {len(df.columns)} cols)")
    print(f"  groups: {sorted(set(sheets['univariate_gwas.csv']['group']))}")


if __name__ == "__main__":
    main()
