#!/usr/bin/env python
"""Prospective temporal-holdout pipeline for the PU immune-target model (Fig 6).

Positives (train): approved for an immune indication by T0 (positives_dated.csv).
Emergent (test): model-target genes whose first interventional trial, any condition,
falls after T0 and that are not in clinic by T0.
Known-non-approved: in clinic by T0 via trial date or approval; held out.

See REPRODUCIBILITY.md, Temporal holdout (T0) labels, for the choice of T0 and the
emergent-gene labelling rule.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find sibling modules
import retrain_at_T0 as rt
import recovery_eval as rv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OT_EVIDENCE = os.path.join(ROOT, "data", "data_drug",
                           "all_drugs_approved_and_in_trial_evidence_ot.csv")
MATRIX = os.path.join(ROOT, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
POSITIVES = os.path.join(ROOT, "data", "prospective", "positives_dated.csv")
DATES = os.path.join(ROOT, "data", "prospective", "drug_trial_dates.json")
# Only tables consumed by the figure script are written to figure_data/; other outputs go to outputs/ as scratch.
FIGDATA = os.path.join(ROOT, "figure_data")
SCRATCH = os.path.join(ROOT, "outputs", "prospective")
os.makedirs(FIGDATA, exist_ok=True)
os.makedirs(SCRATCH, exist_ok=True)

T0_HEADLINE = 2014
T0_SWEEP = [2013, 2014, 2015, 2016, 2017]
# Year-stamped output filenames are derived from T0_HEADLINE, not hardcoded, avoiding stale filenames.
SCORES_CSV = f"scores_at_T0_{T0_HEADLINE}.csv"
EMERGENT_CSV = f"emergent_targets_T0_{T0_HEADLINE}.csv"
SEED = 1234                 # permutation/bootstrap only; the model is seeded by retrain_at_T0


def load_dates():
    """Return (first_trial, approval_year, targets).

    first_trial: per-gene first interventional-trial year across all drugs and conditions.
    approval_year: per-gene immune approval year.
    targets: model-target genes.

    See REPRODUCIBILITY.md, Temporal holdout (T0) labels, for the labelling rule.
    """
    ev = pd.read_csv(OT_EVIDENCE)
    dd = json.load(open(DATES))
    any_year = {c: dd[c].get("any_year") for c in dd}

    # Multi-gene rows in the drug-gene edge table are exploded into one row per gene.
    ev2 = ev[["chembl_id", "gene_target"]].dropna().copy()
    ev2 = ev2.assign(gene=ev2.gene_target.str.split(";")).explode("gene")
    ev2["gene"] = ev2["gene"].str.strip()
    ev2["any_year"] = pd.to_numeric(ev2.chembl_id.map(any_year), errors="coerce")

    first_trial = ev2.groupby("gene")["any_year"].min()

    approval = pd.read_csv(POSITIVES).set_index("gene")["approval_year"].to_dict()

    m = pd.read_parquet(MATRIX)
    targets = set(m.gene[m.pu_role.isin(["P", "trial_heldout"])])
    return first_trial, approval, targets


# ---------------------------------------------------------------------------
# 2. Label sets at a given T0
# ---------------------------------------------------------------------------
def build_split(T0, ft, appr, targets):
    """Partition the labelled target genes at freeze year T0.

    A single first-trial date decides both branches of the split. See REPRODUCIBILITY.md,
    Temporal holdout (T0) labels, for the rule.
    """
    pos = {g for g in targets if pd.notna(appr.get(g)) and appr[g] <= T0}
    emergent, known, dropped = set(), set(), set()
    for g in targets:
        if g in pos:
            continue
        ft_g = ft.get(g, np.nan)
        ap = appr.get(g, np.nan)
        in_clinic = (pd.notna(ft_g) and ft_g <= T0) or (pd.notna(ap) and ap <= T0)
        if in_clinic:
            known.add(g)
        elif pd.notna(ft_g) and ft_g > T0:
            emergent.add(g)
        else:
            # Genes with no trial date or approval are passed as holdout_genes, excluded from negatives and emergent scoring.
            dropped.add(g)
    return pos, emergent, known, dropped


# ---------------------------------------------------------------------------
# 3. Ranking AUC and fixed-model permutation null
# ---------------------------------------------------------------------------
def auc_perm(df, col, n_perm=10000, seed=SEED):
    pool = df[df.in_pool]
    e = pool.is_emergent.values
    s = np.where(pool[col].isna(), np.nanmin(pool[col]) - 1, pool[col].values)
    U, _ = mannwhitneyu(s[e], s[~e], alternative="greater")
    auc = U / (e.sum() * (~e).sum())
    rng = np.random.RandomState(seed)
    order = np.argsort(-s)
    ranks = np.empty_like(order); ranks[order] = np.arange(len(s))
    ne = int(e.sum()); null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.choice(len(s), ne, replace=False)
        null[i] = 1 - (ranks[idx].mean() - (ne - 1) / 2) / (len(s) - ne)
    p = (np.sum(null >= auc) + 1) / (n_perm + 1)
    return auc, p


# ---------------------------------------------------------------------------
# 3b. Paired full-vs-genetics AUC significance on the emergent-vs-unlabeled pool:
#      fast DeLong, paired gene bootstrap, and label-swap permutation.
# ---------------------------------------------------------------------------
def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2


def _fast_delong(preds, m):
    """preds: (k, n) array; the first m columns are the positive-class predictions. Returns (aucs, cov)."""
    k, n = preds.shape; n_neg = n - m
    tx = np.empty((k, m)); ty = np.empty((k, n_neg)); tz = np.empty((k, n))
    for r in range(k):
        tx[r] = _midrank(preds[r, :m]); ty[r] = _midrank(preds[r, m:]); tz[r] = _midrank(preds[r])
    aucs = (tz[:, :m].sum(axis=1) / m - (m + 1) / 2) / n_neg
    v01 = (tz[:, :m] - tx) / n_neg
    v10 = 1 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n_neg
    return aucs, np.atleast_2d(cov)


def delong_paired(scores):
    """Paired DeLong z-test for full-vs-genetics; used by the T0 sweep and headline; no bootstrap/permutation."""
    from scipy.stats import norm
    pool = scores[scores.in_pool].copy()
    y = pool.is_emergent.values.astype(int)
    sf = np.where(pool.score_full.isna(), np.nanmin(pool.score_full) - 1, pool.score_full.values)
    sg = np.where(pool.score_genetics.isna(), np.nanmin(pool.score_genetics) - 1, pool.score_genetics.values)
    n_pool = len(pool); m = int(y.sum())
    order = np.argsort(-y)                      # positives first
    aucs, cov = _fast_delong(np.vstack([sf[order], sg[order]]), m)
    auc_f, auc_g = float(aucs[0]), float(aucs[1])
    delta = auc_f - auc_g
    se = float(np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]))
    z = delta / se
    return dict(n_emergent=m, n_pool=n_pool, auc_full=auc_f, auc_genetics=auc_g,
                delta=delta, se_delta=se, z=z,
                p_two_sided=float(2 * norm.sf(abs(z))), p_one_sided=float(norm.sf(z)),
                ci95_lo=delta - 1.96 * se, ci95_hi=delta + 1.96 * se)


def paired_auc_significance(scores, n_boot=10000, n_perm=10000, seed=SEED):
    pool = scores[scores.in_pool].copy()
    y = pool.is_emergent.values.astype(int)
    sf = np.where(pool.score_full.isna(), np.nanmin(pool.score_full) - 1, pool.score_full.values)
    sg = np.where(pool.score_genetics.isna(), np.nanmin(pool.score_genetics) - 1, pool.score_genetics.values)
    n_pool = len(pool); m = int(y.sum())
    dl = delong_paired(scores)
    auc_f, auc_g = dl["auc_full"], dl["auc_genetics"]
    delta = dl["delta"]; se = dl["se_delta"]; z = dl["z"]
    p_two = dl["p_two_sided"]; p_one = dl["p_one_sided"]

    # paired gene bootstrap: resamples genes and recomputes both AUCs on each draw
    from sklearn.metrics import roc_auc_score
    rng = np.random.RandomState(seed)
    boot = np.empty(n_boot)
    idx_all = np.arange(n_pool)
    for i in range(n_boot):
        bi = rng.choice(idx_all, n_pool, replace=True)
        yb = y[bi]
        if yb.sum() == 0 or yb.sum() == len(yb):
            boot[i] = np.nan; continue
        boot[i] = roc_auc_score(yb, sf[bi]) - roc_auc_score(yb, sg[bi])
    boot = boot[~np.isnan(boot)]
    b_lo, b_hi = np.percentile(boot, [2.5, 97.5])
    p_boot_one = (np.sum(boot <= 0) + 1) / (len(boot) + 1)

    # label-swap permutation: shuffles emergent labels and recomputes the paired AUC delta
    perm = np.empty(n_perm)
    for i in range(n_perm):
        yp = rng.permutation(y)
        perm[i] = roc_auc_score(yp, sf) - roc_auc_score(yp, sg)
    p_perm = (np.sum(np.abs(perm) >= abs(delta)) + 1) / (n_perm + 1)

    return pd.DataFrame([
        {"test": "DeLong (paired correlated ROC)", "n_emergent": m, "n_pool": n_pool,
         "auc_full": round(auc_f, 4), "auc_genetics": round(auc_g, 4), "delta": round(delta, 4),
         "se_delta": round(se, 4), "stat": f"z={z:.3f}", "p_two_sided": round(p_two, 4),
         "p_one_sided": round(p_one, 4), "ci95_lo": round(delta - 1.96 * se, 4),
         "ci95_hi": round(delta + 1.96 * se, 4)},
        {"test": "paired gene bootstrap (10k)", "n_emergent": m, "n_pool": n_pool,
         "auc_full": round(auc_f, 4), "auc_genetics": round(auc_g, 4), "delta": round(float(boot.mean()), 4),
         "se_delta": round(float(boot.std(ddof=1)), 4), "stat": "resample", "p_two_sided": np.nan,
         "p_one_sided": round(float(p_boot_one), 4), "ci95_lo": round(float(b_lo), 4),
         "ci95_hi": round(float(b_hi), 4)},
        {"test": "label-swap permutation (10k)", "n_emergent": m, "n_pool": n_pool,
         "auc_full": round(auc_f, 4), "auc_genetics": round(auc_g, 4), "delta": round(delta, 4),
         "se_delta": np.nan, "stat": "perm", "p_two_sided": round(float(p_perm), 4),
         "p_one_sided": np.nan, "ci95_lo": np.nan, "ci95_hi": np.nan},
    ])


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ft, appr, targets = load_dates()

    # headline split
    pos, emg, known, dropped = build_split(T0_HEADLINE, ft, appr, targets)
    print(f"T0={T0_HEADLINE}: P={len(pos)} emergent={len(emg)} "
          f"known_heldout={len(known)} dropped={len(dropped)}")

    # See REPRODUCIBILITY.md, Temporal holdout (T0) labels.
    full, meta_f = rt.run_T0(pos, emg, T0_HEADLINE, holdout_genes=known | dropped, tag="full")
    gen, meta_g = rt.run_T0(pos, emg, T0_HEADLINE, feature_subset=rt.GENETICS_FEATS,
                            holdout_genes=known | dropped, tag="genetics")
    scores = full.merge(gen[["gene", "score_genetics", "rank_pctile_genetics"]], on="gene")
    scores.to_csv(os.path.join(FIGDATA, SCORES_CSV), index=False)

    # Computes recovery@k with confidence intervals and a null comparison.
    rec = pd.concat([
        rv.recovery_with_ci(scores, "score_full").assign(model="full"),
        rv.recovery_with_ci(scores, "score_genetics").assign(model="genetics"),
    ], ignore_index=True)
    # Named t0_recovery_at_k, not recovery_at_k: panelD_recovery.csv already holds a different recovery analysis (Table 7).
    rec.to_csv(os.path.join(FIGDATA, "t0_recovery_at_k.csv"), index=False)

    # Headline AUC computation.
    af, pf = auc_perm(scores, "score_full")
    ag, pg = auc_perm(scores, "score_genetics")
    auc_df = pd.DataFrame([{"model": "full", "auc": af, "p_emp": pf, "n_emergent": len(emg)},
                           {"model": "genetics", "auc": ag, "p_emp": pg, "n_emergent": len(emg)}])
    auc_df.to_csv(os.path.join(FIGDATA, "auc_emergent.csv"), index=False)
    print(f"AUC full={af:.3f} (p={pf:.4f})  genetics={ag:.3f} (p={pg:.4f})")

    # Paired comparison of full versus genetics feature sets.
    rv.paired_bootstrap_delta(scores, "score_full", "score_genetics").to_csv(
        os.path.join(SCRATCH, "paired_full_vs_genetics.csv"), index=False)

    # Paired AUC significance (DeLong, bootstrap, permutation) is regenerated each run from current scores, so the CSV cannot drift from the figure.
    sig = paired_auc_significance(scores)
    sig.to_csv(os.path.join(SCRATCH, "delong_full_vs_genetics.csv"), index=False)
    d0 = sig.iloc[0]
    print(f"DeLong full-vs-genetics: dAUC={d0.delta:+.4f} {d0.stat} "
          f"p_two={d0.p_two_sided} p_one={d0.p_one_sided}")

    # See REPRODUCIBILITY.md, Temporal holdout permutation resolution.
    rows = []
    for T0 in T0_SWEEP:
        p, e, k, d = build_split(T0, ft, appr, targets)
        # Uses k | d ordering to match the headline analysis, keeping results comparable.
        sc, _ = rt.run_T0(p, e, T0, holdout_genes=k | d, tag="full")
        gsc, _ = rt.run_T0(p, e, T0, feature_subset=rt.GENETICS_FEATS,
                           holdout_genes=k | d, tag="genetics")
        sc = sc.merge(gsc[["gene", "score_genetics", "rank_pctile_genetics"]], on="gene")
        a, pp = auc_perm(sc, "score_full", n_perm=5000)
        dl = delong_paired(sc)
        rows.append(dict(T0=T0, n_P=len(p), n_emergent=len(e), auc_full=a, p_full=pp,
                         auc_genetics=round(dl["auc_genetics"], 4),
                         delta=round(dl["delta"], 4), se_delta=round(dl["se_delta"], 4),
                         delong_p_two_sided=round(dl["p_two_sided"], 4),
                         delong_p_one_sided=round(dl["p_one_sided"], 4),
                         ci95_lo=round(dl["ci95_lo"], 4), ci95_hi=round(dl["ci95_hi"], 4)))
        print(f"  sweep T0={T0}: P={len(p)} E={len(e)} AUC={a:.3f} p={pp:.4f} "
              f"| dAUC={dl['delta']:+.4f} DeLong p2={dl['p_two_sided']:.4f}")
    sens = pd.DataFrame(rows)
    sens.to_csv(os.path.join(FIGDATA, "t0_sensitivity.csv"), index=False)

    # Holds only panel-level summaries (unlabeled count, two medians, worst sweep p); per-gene rank percentiles remain in scores_at_T0_*.csv.
    pool = scores[scores["in_pool"].astype(bool)]
    is_e = pool["is_emergent"].astype(bool)
    e_vals = pool.loc[is_e, "rank_pctile_full"].dropna()
    u_vals = pool.loc[~is_e, "rank_pctile_full"].dropna()
    eg_vals = pool.loc[is_e, "rank_pctile_genetics"].dropna()
    # Paired subset requires both rank percentiles present per gene; e_vals/eg_vals drop NaNs independently, so their n can exceed the pairable set.
    paired = pool.loc[is_e, ["rank_pctile_full", "rank_pctile_genetics"]].dropna()
    pd.DataFrame([{
        "median_emergent": float(np.median(e_vals)),
        "median_unlabeled": float(np.median(u_vals)),
        "n_unlabeled": int(len(u_vals)),
        "p_full_max": float(sens["p_full"].max()),
        "median_emergent_genetics": float(np.median(eg_vals)),
        "q1_emergent": float(np.quantile(e_vals, 0.25)),
        "q3_emergent": float(np.quantile(e_vals, 0.75)),
        "q1_emergent_genetics": float(np.quantile(eg_vals, 0.25)),
        "q3_emergent_genetics": float(np.quantile(eg_vals, 0.75)),
        "n_emergent_paired": int(len(paired)),
        "n_full_higher": int((paired["rank_pctile_full"] > paired["rank_pctile_genetics"]).sum()),
    }]).to_csv(os.path.join(FIGDATA, "t0_panel_stats.csv"), index=False)

    # Emergent target table.
    et = pd.DataFrame({"gene": sorted(emg)})
    # Date column records the split date that assigned the gene to the test set.
    et["first_trial_year"] = et.gene.map(ft).astype(int)
    et = et.merge(scores[["gene", "rank_pctile_full", "rank_pctile_genetics"]], on="gene")
    et = et.sort_values("rank_pctile_full", ascending=False)
    et.to_csv(os.path.join(FIGDATA, EMERGENT_CSV), index=False)

    # Split JSON.
    split = dict(
        T0=T0_HEADLINE,
        n_pos_at_T0=len(pos), n_emergent=len(emg), n_known_nonapproved_at_T0=len(known),
        n_feat_full=meta_f["n_feat"], n_feat_genetics=meta_g["n_feat"],
        pos_at_T0=sorted(pos), emergent=sorted(emg),
        known_nonapproved_at_T0=sorted(known), dropped_undated=sorted(dropped),
        sweep={int(r.T0): dict(P=int(r.n_P), E=int(r.n_emergent),
                               auc=round(float(r.auc_full), 3), p=float(r.p_full))
               for r in sens.itertuples()},
        design_note=("P=approved-by-T0 (positives_dated.csv); ONE date rule throughout -- the "
                     "CT.gov first-trial year over ALL drugs, any condition, decides both "
                     "sides of the split, so no hand-set threshold enters the design. "
                     "emergent=first-any-trial>T0, model targets, not in-clinic by T0, NO "
                     "immune-condition filter; known_nonapproved=broad fail-safe (any-trial "
                     "or approval <=T0) held out of negatives; dropped=no datable trial and "
                     "no approval, ALSO held out of negatives (both are labelled targets that "
                     "cannot be placed on the timeline, so neither is evidence of "
                     "non-targetness; matches the production model, which never draws a "
                     "trial_heldout gene as a pseudo-negative). Adopted 2026-08-28, replacing "
                     "a design that dated emergent "
                     "membership from selective drugs only; consequence: 78 ribosomal proteins "
                     "dated solely via MT-3724 (anti-CD20/Shiga-toxin conjugate) enter the "
                     "T0<=2014 test set; being strongly constrained, they favour the "
                     "genetics-only comparator and shrink the FG increment."),
    )
    json.dump(split, open(os.path.join(FIGDATA, "t0_split.json"), "w"), indent=1)

    print("wrote data tables to", FIGDATA, "and", SCRATCH)
    print("render the figure with:  Rscript src/figures/fig_temporal_holdout.R")


if __name__ == "__main__":
    main()
