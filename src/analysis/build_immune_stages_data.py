#!/usr/bin/env python3
"""Stage-1 tables for the immune-indication trial-validation figure (Fig 5), rendered by
src/figures/fig_immune_stages.R.

Panel a: top-50 non-approved PU nominations, barred by the highest clinical stage of an
immune-indication drug against that target.
Panel b: fraction of genes with an immune-indication drug among the full model's top 50 and
the genetics-only model's top 50, each against a bootstrapped control, for any immune trial
(Ph I+) and advanced trials (Ph III), with 95% bootstrap intervals and enrichment p-values.
Panel c: the any-immune-trial arm of panel b at three shortlist depths (25/50/100).

See REPRODUCIBILITY.md, Other significance tests, for the statistical caveats on both panels.

Deterministic: SEED=1234, N_TOP=N_CTRL=50, B=10000.

Inputs (committed; no live API calls):
  - outputs/model/full_model_pu_scores.csv               PU score and role per gene
  - outputs/model/genetics_only_pu_scores.csv            genetics-only score and role per gene
  - data/data_drug/all_drugs_approved_and_in_trial_evidence_ot.csv, column
    drug_stage_for_immune_indication
  - figure_data/ranked_atlas.csv                          reporting columns merged onto top 50

Outputs:
  - figure_data/immune_stages_top50.csv     top 50 with reporting columns
  - figure_data/immune_stages_counts.json   observed/control rates, CIs, folds, p-values

Runs after build_ranked_atlas_table.py, which writes the atlas this script reads.
"""
import os
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _stats import bh  # noqa: E402

# ---------------------------------------------------------------- config
SEED, N_TOP, N_CTRL, B, CTRL_RANK_MIN = 1234, 50, 50, 10000, 200
# The control draw size always equals the shortlist depth; a mismatched size would misstate the interval width and could flatter the comparison.
DEPTHS_C = (25, N_TOP, 100)
assert max(DEPTHS_C) <= CTRL_RANK_MIN, (
    f"rank floor {CTRL_RANK_MIN} must exclude the whole of every panel-c depth {DEPTHS_C}: above "
    "the floor the control pool would contain nominated genes and the null would stop being a "
    "background")
assert N_TOP in DEPTHS_C, "panel c must include panel b's depth or the two panels stop comparing"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
FIGDATA = os.path.join(ROOT, "figure_data")
os.makedirs(FIGDATA, exist_ok=True)
PU_SCORES = os.path.join(ROOT, "outputs", "model", "full_model_pu_scores.csv")
GEN_SCORES = os.path.join(ROOT, "outputs", "model", "genetics_only_pu_scores.csv")
EVIDENCE  = os.path.join(ROOT, "data", "data_drug", "all_drugs_approved_and_in_trial_evidence_ot.csv")
ATLAS     = os.path.join(FIGDATA, "ranked_atlas.csv")
OUT_CSV   = os.path.join(FIGDATA, "immune_stages_top50.csv")
OUT_JSON  = os.path.join(FIGDATA, "immune_stages_counts.json")

# Threshold matches pu_target_model.py's top-novel table; unmeasured genes are not_measured, not weak.
DONOR_VETTED_MIN = 0.3

# STAGE_NUM in fig_immune_stages.R must match this label encoding; a mismatch shifts bar lengths silently.
STAGE_ORD = {"Approved": 4, "Phase 3": 3, "Phase 2": 2, "Phase 1": 1, "Unknown": 0.5}
NUM2LAB   = {4: "Approved", 3: "Phase 3", 2: "Phase 2", 1: "Phase 1", 0.5: "Unknown"}

# ---------------------------------------------------------------- data
g = pd.read_csv(PU_SCORES)
gen = pd.read_csv(GEN_SCORES)
ev = pd.read_csv(EVIDENCE)

# Split gene_target on ';' and explode before grouping; else ~2,617/7,030 rows undercount via composite keys.
ev["stg_num"] = ev["drug_stage_for_immune_indication"].map(STAGE_ORD)
ev_g = ev.assign(gene_target=ev["gene_target"].str.split(";")).explode("gene_target")
imm_num_by_gene = ev_g.groupby("gene_target")["stg_num"].max()
imm_lab_by_gene = imm_num_by_gene.map(NUM2LAB)

# >=1 excludes the Unknown stage (0.5): a drug with no recorded stage is not evidence of a trial.
phI  = lambda df: df.imm_num >= 1          # any immune trial (Ph I+)
ph3  = lambda df: df.imm_num >= 3          # advanced immune trial (Ph III only)

ci   = lambda a: (float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))
pval = lambda rates, hits, n_top: float((rates >= 100 * hits / n_top).mean())
fold = lambda rates, hits, n_top: float((hits / n_top) / (rates.mean() / 100))


def rank_nonapproved(scores, score_col):
    """Non-approved pool (trial + unlabeled) ranked by score_col, with immune stage attached. pu_role is shared across both score files, so the pool is the same 19,162 genes for both models, only the order differs."""
    na = (scores[scores.pu_role != "P"]
          .sort_values(score_col, ascending=False).reset_index(drop=True))
    na["rank_na"] = np.arange(1, len(na) + 1)
    na["furthest_imm_stage"] = na.gene.map(imm_lab_by_gene)
    na["imm_num"] = na.gene.map(imm_num_by_gene).fillna(-1)
    return na


def bootstrap_arm(na, n_top):
    """Observed top-n rates and bootstrapped control for one model at one depth. See REPRODUCIBILITY.md, Other significance tests, for the control-sampling rationale."""
    top_ = na.head(n_top)
    t_a, t_d = int(phI(top_).sum()), int(ph3(top_).sum())
    pool_ = na[na.rank_na > CTRL_RANK_MIN]
    p_a_v, p_d_v = phI(pool_).values, ph3(pool_).values
    rng_, idx_ = np.random.RandomState(SEED), np.arange(len(pool_))
    a_rates, d_rates = np.empty(B), np.empty(B)
    for b in range(B):
        sel = rng_.choice(idx_, size=n_top, replace=True)
        a_rates[b] = p_a_v[sel].mean() * 100
        d_rates[b] = p_d_v[sel].mean() * 100
    return top_, pool_, t_a, t_d, a_rates, d_rates


def arm_stats(na, n_top, prefix):
    """One (model, depth) cell as JSON keys under prefix; panel b's two cells only."""
    top_, pool_, t_a, t_d, a_rates, d_rates = bootstrap_arm(na, n_top)
    return top_, {
        f"{prefix}ctrl_pool_n": int(len(pool_)),
        f"{prefix}top_phI_plus": t_a, f"{prefix}top_ph3_plus": t_d,
        f"{prefix}top_phI_pct": 100 * t_a / n_top, f"{prefix}top_ph3_pct": 100 * t_d / n_top,
        f"{prefix}ctrl_phI_mean_pct": float(a_rates.mean()),
        f"{prefix}ctrl_phI_ci": list(ci(a_rates)),
        f"{prefix}ctrl_ph3_mean_pct": float(d_rates.mean()),
        f"{prefix}ctrl_ph3_ci": list(ci(d_rates)),
        f"{prefix}fold_phI": fold(a_rates, t_a, n_top),
        f"{prefix}fold_ph3": fold(d_rates, t_d, n_top),
        f"{prefix}bootstrap_p_phI": pval(a_rates, t_a, n_top),
        f"{prefix}bootstrap_p_ph3": pval(d_rates, t_d, n_top),
        f"{prefix}top_immune_genes": top_[top_.imm_num > 0][["gene", "furthest_imm_stage"]].to_dict("records"),
    }


def sweep_record(na, n_top, model):
    """One (model, depth) row of depth_sweep, flat and scalar-valued for the renderer. CI bounds are separate lo/hi fields since jsonlite cannot read a nested list column into a data.frame."""
    _, pool_, t_a, t_d, a_rates, d_rates = bootstrap_arm(na, n_top)
    lo_a, hi_a = ci(a_rates)
    lo_d, hi_d = ci(d_rates)
    return {
        "depth": n_top, "model": model, "ctrl_pool_n": int(len(pool_)),
        "top_phI_plus": t_a, "top_phI_pct": 100 * t_a / n_top,
        "ctrl_phI_mean_pct": float(a_rates.mean()),
        "ctrl_phI_ci_lo": lo_a, "ctrl_phI_ci_hi": hi_a,
        "fold_phI": fold(a_rates, t_a, n_top), "bootstrap_p_phI": pval(a_rates, t_a, n_top),
        # Written but not plotted in panel c; see the panel c note in the module docstring.
        "top_ph3_plus": t_d, "top_ph3_pct": 100 * t_d / n_top,
        "ctrl_ph3_mean_pct": float(d_rates.mean()),
        "ctrl_ph3_ci_lo": lo_d, "ctrl_ph3_ci_hi": hi_d,
        "fold_ph3": fold(d_rates, t_d, n_top), "bootstrap_p_ph3": pval(d_rates, t_d, n_top),
    }


# ---------------------------------------------------------------- stats
# Panel b: both thresholds at depth 50, full model and genetics-only.
nonappr = rank_nonapproved(g, "pu_score")
nonappr_gen = rank_nonapproved(gen, "genetics_score")

top, stats_b_full = arm_stats(nonappr, N_TOP, "")
top_gen, stats_b_gen = arm_stats(nonappr_gen, N_TOP, "gen_")

# Panel c: the same two arms across DEPTHS_C.
sweep = [sweep_record(na, d, m)
         for d in DEPTHS_C
         for na, m in ((nonappr, "full"), (nonappr_gen, "genetics_only"))]

# Depth 50 is shared by panels b and c; bootstrap_arm re-seeds per call so both draw identical samples, and this assertion catches any drift.
for rec in sweep:
    if rec["depth"] != N_TOP:
        continue
    src = stats_b_full if rec["model"] == "full" else stats_b_gen
    pre = "" if rec["model"] == "full" else "gen_"
    for f_ in ("top_phI_plus", "top_phI_pct", "ctrl_phI_mean_pct", "fold_phI", "bootstrap_p_phI",
               "top_ph3_plus", "fold_ph3", "bootstrap_p_ph3"):
        assert rec[f_] == src[pre + f_], \
            f"panel b and panel c disagree on {rec['model']} {f_} at depth {N_TOP}"

# ---------------------------------------------------------------- multiple-testing correction
# BH applied across the eight enrichment tests Figure 5 displays, as one family. See REPRODUCIBILITY.md, Multiple-testing correction, for the rationale.
sites = [(stats_b_full, "bootstrap_p_phI"), (stats_b_full, "bootstrap_p_ph3"),
         (stats_b_gen, "gen_bootstrap_p_phI"), (stats_b_gen, "gen_bootstrap_p_ph3")]
sites += [(rec, "bootstrap_p_phI") for rec in sweep if rec["depth"] != N_TOP]
assert len(sites) == 8, f"Figure 5 displays eight enrichment tests, collected {len(sites)}"
for (container, key), q in zip(sites, bh([c[k] for c, k in sites])):
    container[key + "_bh"] = q

# Depth-50 sweep rows reuse panel b's adjusted value rather than being corrected again, so every sweep record shares one key for jsonlite.
for rec in sweep:
    if rec["depth"] == N_TOP:
        src = stats_b_full if rec["model"] == "full" else stats_b_gen
        pre = "" if rec["model"] == "full" else "gen_"
        rec["bootstrap_p_phI_bh"] = src[pre + "bootstrap_p_phI_bh"]

# The two model bars in a panel are paired by matched genes at each depth, not independent samples; their difference is not a two-sample test.
shared = [{"depth": d,
           "n_shared": int(len(set(nonappr.gene.head(d)) & set(nonappr_gen.gene.head(d))))}
          for d in DEPTHS_C]
n_shared = next(r["n_shared"] for r in shared if r["depth"] == N_TOP)

# ---------------------------------------------------------------- data out
# Reporting columns joined from the atlas via how="left" on a validated one-to-one key, so a missing gene appears as NaN rather than shortening the table; checked by the assertion below.
atlas = pd.read_csv(ATLAS, usecols=["gene", "pu_score", "rank", "rank_pctile",
                                    "crossdonor_confidence"])
sheet = top[["rank_na", "gene", "pu_role", "furthest_imm_stage"]].rename(
    columns={"rank_na": "rank_nonapproved"})
sheet = sheet.merge(atlas, on="gene", how="left", validate="one_to_one")
assert sheet.notna()["rank"].all(), \
    f"not in ranked_atlas.csv: {sorted(sheet.loc[sheet['rank'].isna(), 'gene'])}"

# See PROVENANCE.md, Trial-validation figure — `final_plots/figure_trial_validation_immune_stages_R.png` for the bit-identical scores check.
assert (sheet["pu_score"].values == top["pu_score"].values).all(), \
    "pu_score disagrees between ranked_atlas.csv and full_model_pu_scores.csv"

# The stage ships once, as furthest_imm_stage, recomputed here from the drug evidence table. The
# ranked atlas carries the same annotation via pu_model_matrix.parquet; it used to be reported
# beside this one and asserted equal, which was the only check that the two sources had not drifted.

cd = sheet["crossdonor_confidence"]
sheet["donor_vetted"] = np.where(cd.notna(),
                                 np.where(cd >= DONOR_VETTED_MIN, "yes", "weak"), "not_measured")

sheet = sheet[["rank_nonapproved", "gene", "pu_score", "rank", "rank_pctile", "pu_role",
               "furthest_imm_stage", "crossdonor_confidence", "donor_vetted"]]
assert len(sheet) == N_TOP, f"expected {N_TOP} nominations, got {len(sheet)}"
sheet.to_csv(OUT_CSV, index=False)
counts = {
    "source": "curated Open Targets: all_drugs_approved_and_in_trial_evidence_ot.csv (immune stage per drug)",
    "immune_stage_metric": "per-gene furthest stage of any immune-indication drug",
    "N_top": N_TOP, "N_ctrl_per_draw": N_CTRL, "n_bootstrap": B,
    "ctrl_pool_rank_min": CTRL_RANK_MIN, "ctrl_pool_n": stats_b_full["ctrl_pool_n"], "seed": SEED,
    "panelC_depths": list(DEPTHS_C),
    # Unprefixed keys are the full model at depth 50 (panel b); gen_ is the genetics-only arm. *_top_immune_genes use >0, retaining Unknown (0.5); phI/ph3 counts use >=1 and exclude Unknown.
    **{k: v for k, v in stats_b_full.items() if k != "ctrl_pool_n"}, **stats_b_gen,
    # Panel c plots only the phI columns of these rows; ph3 columns are written but not plotted.
    "depth_sweep": sweep,
    "shared_by_depth": shared,
    "n_top50_shared": n_shared,
}
json.dump(counts, open(OUT_JSON, "w"), indent=2)

for rec in sweep:
    star = "  <- panel b" if rec["depth"] == N_TOP else ""
    print(f"top{rec['depth']:3d} {rec['model']:13s}: "
          f"PhI+={rec['top_phI_plus']:2d} ({rec['top_phI_pct']:.0f}%) "
          f"fold {rec['fold_phI']:.1f}x p={rec['bootstrap_p_phI']:.4f} | ctrl "
          f"{rec['ctrl_phI_mean_pct']:.2f}%  [not plotted: Ph3={rec['top_ph3_plus']} "
          f"fold {rec['fold_ph3']:.1f}x p={rec['bootstrap_p_ph3']:.4f}]{star}")
print("shared genes by depth: " + ", ".join(f"top{r['depth']}={r['n_shared']}" for r in shared))
print(f"wrote {OUT_CSV} and {OUT_JSON}")
