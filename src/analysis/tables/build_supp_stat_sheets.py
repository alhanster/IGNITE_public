#!/usr/bin/env python3
"""Flattens statistics JSON outputs into supplementary-table sheets under figure_data/, covering Supp Tables 7, 8, 10, 11 and 14.

See REPRODUCIBILITY.md for the step order.
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import _platform  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"

FIGDATA = os.path.join(ROOT, "figure_data")


def fd(name):
    return os.path.join(FIGDATA, name)


def load_json(name):
    with open(fd(name)) as fh:
        return json.load(fh)


# --------------------------------------------------------------------------------------------
# Supp Table 7: discrimination statistics
# --------------------------------------------------------------------------------------------
def heldout_trial_stats(j):
    """heldout_trial_auc.json becomes a long tidy frame since the statistics are heterogeneous; value stays object dtype so counts like n=727 are not rendered as 727.0."""
    rows = [
        ("n_in_trial", "", j["n_in_trial"], None, None),
        ("n_non_target", "", j["n_non_target"], None, None),
        ("n_features", "full", j["n_feat_full"], None, None),
        ("n_features", "genetics", j["n_feat_genetics"], None, None),
        ("auc", "full", j["auc_full"], *j["auc_full_ci"]),
        ("auc", "genetics", j["auc_genetics"], *j["auc_genetics_ci"]),
        ("auc_delta", "full - genetics", j["delta"], *j["delta_ci95"]),
        ("delong_z", "full vs genetics", j["delong_z"], None, None),
        ("delong_p", "full vs genetics", j["delong_p"], None, None),
    ]
    # dtype=object is set at construction, not via astype() afterward, since by then pandas has already widened the column to float64 and the integers are gone.
    return pd.DataFrame(rows, columns=["statistic", "model", "value", "ci_lo", "ci_hi"],
                        dtype=object)


# Group order and naming: see REPRODUCIBILITY.md.
GROUP_ORDER = ["approved", "in-trial", "non-target"]


def score_by_group_summary(j):
    g = j["groups"]
    return pd.DataFrame([{"group": k, **{s: g[k][s] for s in
                                         ("n", "median", "mean", "q25", "q75")}}
                         for k in GROUP_ORDER])


def score_by_group_tests(j):
    """Omnibus tests, the headline held-out test, and the three pairwise comparisons; stars mixes string significance annotations into a numeric column so each row matches the annotation on its figure panel."""
    ho = j["heldout_in_trial_vs_non_target"]
    rows = [
        {"test": "Kruskal-Wallis", "group_a": "", "group_b": "", "statistic": "H",
         "value": j["kruskal_H"], "p": j["kruskal_p"], "p_holm": None, "stars": ""},
        {"test": "one-way ANOVA", "group_a": "", "group_b": "", "statistic": "F",
         "value": j["anova_F"], "p": j["anova_p"], "p_holm": None, "stars": ""},
        {"test": "Mann-Whitney (one-sided)", "group_a": "in-trial", "group_b": "non-target",
         "statistic": "U", "value": ho["U"], "p": ho["p_one_sided"], "p_holm": None, "stars": ""},
    ]
    for pw in j["pairwise_mannwhitney"]:
        rows.append({"test": "Mann-Whitney (two-sided)", "group_a": pw["a"], "group_b": pw["b"],
                     "statistic": "U", "value": pw["U"], "p": pw["p"],
                     "p_holm": pw["p_holm"], "stars": pw["stars"]})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------
# Supp Table 8: scrambled-feature control
# --------------------------------------------------------------------------------------------
def scrambled_control_stats(j):
    """One long frame of the observed increment, the 99 draws, then summary scalars. key is the draw index on scrambled_draw rows and a statistic name on summary rows; draws sort numerically, since the JSON keys are strings and would otherwise order 0, 1, 10, 11. value stays object dtype so counts like n_scramble and n_ge_real are not rendered as floats."""
    rows = [{"kind": "observed", "key": "real_increment", "value": j["real"]}]
    rows += [{"kind": "scrambled_draw", "key": str(i), "value": j["scrambled"][str(i)]}
             for i in sorted(int(k) for k in j["scrambled"])]
    s = j["summary"]
    rows += [{"kind": "summary", "key": k, "value": s[k]} for k in
             ("n_scramble", "scrambled_mean", "scrambled_sd", "scrambled_max", "n_ge_real",
              "emp_p")]
    return pd.DataFrame(rows, dtype=object)


# --------------------------------------------------------------------------------------------
# Supp Table 10: bootstrap trial enrichment
# --------------------------------------------------------------------------------------------
def immune_stages_bootstrap(j):
    """See REPRODUCIBILITY.md."""
    rows = []
    for model, pre in [("full", ""), ("genetics_only", "gen_")]:
        for label, key in [("Phase >= 1 (any immune trial)", "phI"),
                           ("Phase >= 3 (advanced immune trial)", "ph3")]:
            obs_pct, ctrl_pct = j[f"{pre}top_{key}_pct"], j[f"{pre}ctrl_{key}_mean_pct"]
            lo, hi = j[f"{pre}ctrl_{key}_ci"]
            rows.append({
                "model": model, "threshold": label,
                "observed_count": j[f"{pre}top_{key}_plus"], "observed_pct": obs_pct,
                "expected_count": ctrl_pct * j["N_top"] / 100, "expected_pct": ctrl_pct,
                "ci95_lo_pct": lo, "ci95_hi_pct": hi,
                "fold_enrichment": j[f"{pre}fold_{key}"],
                "p_bootstrap": j[f"{pre}bootstrap_p_{key}"],
                "p_bootstrap_bh": j[f"{pre}bootstrap_p_{key}_bh"],
                "n_top": j["N_top"], "n_ctrl_per_draw": j["N_ctrl_per_draw"],
                "n_bootstrap": j["n_bootstrap"], "ctrl_pool_rank_min": j["ctrl_pool_rank_min"],
                "ctrl_pool_n": j[f"{pre}ctrl_pool_n"], "seed": j["seed"]})
    # dtype=object because the six trailing count columns would otherwise render as floats, e.g. 50.0.
    return pd.DataFrame(rows, dtype=object)


def immune_stages_depth_sweep(j):
    """See REPRODUCIBILITY.md."""
    rows = []
    for r in j["depth_sweep"]:
        ctrl_pct = r["ctrl_phI_mean_pct"]
        rows.append({
            "depth": r["depth"], "model": r["model"],
            "threshold": "Phase >= 1 (any immune trial)",
            "observed_count": r["top_phI_plus"], "observed_pct": r["top_phI_pct"],
            # round(...,10) removes float noise only: at depth 100, x100/100 does not round-trip; source rates have 4 decimals.
            "expected_count": round(ctrl_pct * r["depth"] / 100, 10),
            "expected_pct": ctrl_pct,
            "ci95_lo_pct": r["ctrl_phI_ci_lo"], "ci95_hi_pct": r["ctrl_phI_ci_hi"],
            "fold_enrichment": r["fold_phI"], "p_bootstrap": r["bootstrap_p_phI"],
            # Depth-50 rows carry panel b's adjusted value, correcting each bar once; only phI has a p_bootstrap_bh entry, since only phI is plotted there.
            "p_bootstrap_bh": r["bootstrap_p_phI_bh"],
            "n_bootstrap": j["n_bootstrap"], "ctrl_pool_rank_min": j["ctrl_pool_rank_min"],
            "ctrl_pool_n": r["ctrl_pool_n"], "seed": j["seed"]})
    # dtype=object here for the same count-column reason as the preceding sheet.
    return pd.DataFrame(rows, dtype=object)


# --------------------------------------------------------------------------------------------
# Supp Table 11: leakage-controlled paired DeLong
# --------------------------------------------------------------------------------------------
def leakage_controlled_delong_stats(j):
    """arm is a column, as in heldout_trial_stats.csv; note (provenance) is dropped, see PROVENANCE.md."""
    rows = []
    for arm in ("standard_taskA", "leakage_controlled"):
        a = j["arms"][arm]
        rows += [(arm, "n_genes", a["n_genes"]),
                 (arm, "n_pos", a["n_pos"]),
                 (arm, "auc_pu_full_model", a["auc_pu"]),
                 (arm, "auc_gps_overall", a["auc_gps"]),
                 (arm, "auc_delta", a["delta"]),
                 (arm, "delong_z", a["z"]),
                 (arm, "delong_p", a["p"]),
                 (arm, "delong_p_holm", a["p_holm"]),
                 (arm, "significance", a["stars"])]
    # n_features_full_model is model-level, not per-arm, so it carries no arm label.
    rows.append(("", "n_features_full_model", j["n_feat_full"]))
    # dtype set at column construction, not via astype() after, or counts like 19162 render as floats.
    return pd.DataFrame(rows, columns=["arm", "statistic", "value"], dtype=object)


# --------------------------------------------------------------------------------------------
# Supp Table 14: STAT4 knockdown-neighborhood summary
# --------------------------------------------------------------------------------------------
def vignette_knn_summary(j):
    """vignette_knn_stats.json: condition and focal record the one condition and anchor gene each row's cosines use."""
    rows = [("focal_gene", j["focal"]),
            ("condition", j["condition"]),
            ("n_anchor_rows", j["n_anchor_rows"]),
            ("n_unlabeled_rows", j["n_unlabeled_rows"]),
            ("knn_p95", j["knn_p95"]),
            ("knn_median_unlabeled", j["knn_median_unlabeled"])]
    return pd.DataFrame(rows, columns=["statistic", "value"], dtype=object)


# --------------------------------------------------------------------------------------------
def check_heldout(j):
    """heldout_trial_auc.json's verification block, tested against the file it names."""
    assert j["n_in_trial"] == 727, f"expected 727 held-out in-trial genes, got {j['n_in_trial']}"

    rec = pd.read_csv(fd("panelD_recovery.csv")).set_index("k")["genetics"]
    claimed = {int(k): int(v) for k, v in j["verification"]["recovery_counts"].items()}
    actual = {int(k): int(rec[k]) for k in claimed}
    assert claimed == actual, (
        "heldout_trial_auc.json claims the genetics-only ranking reproduces panelD_recovery.csv, "
        f"but its recovery_counts {claimed} do not match the file's {actual}")
    print(f"  recovery cross-check: {claimed} matches panelD_recovery.csv")


def check_score_by_group(j):
    """Group sizes in the JSON must be the row counts of the CSV that sits beside it."""
    counts = pd.read_csv(fd("score_by_group.csv"))["grp"].value_counts().to_dict()
    for g in GROUP_ORDER:
        assert j["groups"][g]["n"] == counts.get(g), (
            f"score_by_group_stats.json says n={j['groups'][g]['n']} for {g!r}, but "
            f"score_by_group.csv has {counts.get(g)} rows")
    total = sum(j["groups"][g]["n"] for g in GROUP_ORDER)
    assert total == 19502, f"group sizes sum to {total}, not the 19,502-gene matrix"

    # The two Mann-Whitney rows over the same pair are complementary, not identical; a mismatch means one was computed over a different pool.
    n_it, n_nt = j["groups"]["in-trial"]["n"], j["groups"]["non-target"]["n"]
    pair = next(p for p in j["pairwise_mannwhitney"]
                if {p["a"], p["b"]} == {"in-trial", "non-target"})
    got = j["heldout_in_trial_vs_non_target"]["U"] + pair["U"]
    assert got == n_it * n_nt, f"the two U statistics sum to {got}, not {n_it} x {n_nt}"
    print(f"  group sizes reconcile ({total} genes); the two U statistics sum to {n_it}x{n_nt}")


def check_permutation(j_scr):
    """See REPRODUCIBILITY.md."""
    null = pd.read_csv(fd("label_permutation_null.csv"))
    pv = pd.read_csv(fd("label_permutation_pvalues.csv"))
    ladder = pd.read_csv(fd("panelA_auc_ladder.csv")).set_index("model")["auc"]

    n = len(null)
    assert (pv["n_perm"] == n).all(), f"n_perm disagrees with the {n} rows of the null"
    for _, r in pv.iterrows():
        k = int((null[r["metric"]] >= r["observed"]).sum())
        expect = (k + 1) / (n + 1)
        assert abs(expect - r["emp_p"]) < 1e-12, (
            f"{r['metric']}: emp_p {r['emp_p']} is not the continuity-corrected "
            f"({k}+1)/({n}+1) = {expect}")

    # The four level rows are the ladder's AUCs; equality is exact since both files are written from the same in-memory run.
    for metric, model in [("genetic", "genetic only"),
                          ("genetic+observational", "+ observational"),
                          ("genetic+perturbational", "+ perturbational"), ("full", "full")]:
        obs = float(pv.loc[pv["metric"] == metric, "observed"].iloc[0])
        # The permutation tables come from the opt-in permutation-null target on the reference
        # platform; a refit elsewhere moves the ladder, so there this is reported, not fatal.
        _platform.check(obs == ladder[model], (
            f"{metric}: label_permutation_pvalues.csv says {obs}, panelA_auc_ladder.csv "
            f"says {ladder[model]} -- the two tables describe different runs"))

    # The full-genetic increment appears in three files; float arithmetic gives ~4e-17 discrepancy, so this is a tolerance check.
    # Taken from the two level rows: the increment is no longer a row of its own, the label-permutation null being no test of one rung against another.
    inc = (float(pv.loc[pv["metric"] == "full", "observed"].iloc[0])
           - float(pv.loc[pv["metric"] == "genetic", "observed"].iloc[0]))
    for label, other in [("scrambled_feature_control.json", j_scr["real"]),
                         ("panelA_auc_ladder.csv", ladder["full"] - ladder["genetic only"])]:
        _platform.check(abs(inc - other) < 1e-15,
                        f"full-genetic increment is {inc} here and {other} in {label}")
    print(f"  all {len(pv)} emp_p re-derive from the {n} draws; 4 levels match the ladder exactly")


def check_scrambled(j):
    n = j["summary"]["n_scramble"]
    assert n == 99, f"the published control is 99 draws, this JSON says {n}"
    assert len(j["scrambled"]) == n, \
        f"summary claims {n} draws but the file carries {len(j['scrambled'])}"


def check_immune_stages(j):
    """See REPRODUCIBILITY.md."""
    assert j["N_top"] == 50, f"the published table is a top 50, this JSON says {j['N_top']}"
    top = pd.read_csv(fd("immune_stages_top50.csv"))
    assert len(top) == j["N_top"], f"immune_stages_top50.csv has {len(top)} rows, not {j['N_top']}"

    # Use stages, not notna(): Phase>=1 excludes Unknown-stage genes; they align only if none is Unknown-only.
    staged = int(top["furthest_imm_stage"].isin(["Phase 1", "Phase 2", "Phase 3", "Approved"]).sum())
    assert staged == j["top_phI_plus"], (
        f"the JSON counts {j['top_phI_plus']} nominations at Phase >= 1, but "
        f"immune_stages_top50.csv carries {staged}")
    advanced = int(top["furthest_imm_stage"].isin(["Phase 3", "Approved"]).sum())
    assert advanced == j["top_ph3_plus"], (
        f"the JSON counts {j['top_ph3_plus']} at Phase >= 3, the sheet carries {advanced}")

    # The 9-held-out claim is a sheet property, not the JSON's, and the line reviewers will most likely check.
    roles = set(top.loc[top["furthest_imm_stage"].notna(), "pu_role"])
    assert roles == {"trial_heldout"}, \
        f"nominations with an immune-indication drug should all be trial_heldout, got {roles}"

    for key in ("phI", "ph3"):
        got = (j[f"top_{key}_pct"] / j[f"ctrl_{key}_mean_pct"])
        assert abs(got - j[f"fold_{key}"]) < 1e-12, \
            f"fold_{key} is {j[f'fold_{key}']} but observed/expected is {got}"

    # Unprefixed keys: full-model depth-50 values; gen_ keys: genetics-only. Keeps Table 10 matched to its bootstrap sheet.
    for rec in j["depth_sweep"]:
        if rec["depth"] != j["N_top"]:
            continue
        pre = "" if rec["model"] == "full" else "gen_"
        # Only bootstrap_p_phI_bh exists on sweep records; adding bootstrap_p_ph3_bh here would KeyError.
        for f_ in ("top_phI_plus", "top_phI_pct", "ctrl_phI_mean_pct", "fold_phI",
                   "bootstrap_p_phI", "bootstrap_p_phI_bh",
                   "top_ph3_plus", "top_ph3_pct", "ctrl_ph3_mean_pct",
                   "fold_ph3", "bootstrap_p_ph3"):
            assert rec[f_] == j[pre + f_], (
                f"depth_sweep's depth-{j['N_top']} {rec['model']} row has {f_} = {rec[f_]}, but "
                f"the panel-b key {pre + f_} is {j[pre + f_]} -- Table 10's two bootstrap sheets "
                "would disagree about the same bar")
    n_sweep = len(j["depth_sweep"])
    assert n_sweep == 2 * len(j["panelC_depths"]), \
        f"depth_sweep should be 2 models x {len(j['panelC_depths'])} depths, got {n_sweep} rows"

    # gen_ keys alone hold the genetics arm's depth-50 values; a missing one silently duplicates the full model's.
    for key in ("phI", "ph3"):
        for f_ in (f"top_{key}_plus", f"top_{key}_pct", f"ctrl_{key}_mean_pct", f"ctrl_{key}_ci",
                   f"fold_{key}", f"bootstrap_p_{key}", f"bootstrap_p_{key}_bh"):
            assert f"gen_{f_}" in j, f"immune_stages_counts.json has no gen_{f_}"
            assert j[f"gen_{f_}"] != j[f_] or f_.endswith("_ci"), (
                f"gen_{f_} and {f_} are both {j[f_]} -- the two arms of Fig 5b differ at every "
                "other quantity, so identical values here mean the genetics rows read the full "
                "model's keys")
    assert j["gen_ctrl_pool_n"] == j["ctrl_pool_n"], (
        f"the two control pools are {j['ctrl_pool_n']} and {j['gen_ctrl_pool_n']} genes. They "
        "have always been the same 19,162 non-approved genes in a different order; if that has "
        "changed, Fig 5b's single grey bar and this sheet's per-arm ctrl_pool_n both need review")
    print(f"  bootstrap: {staged}/{j['N_top']} at Phase>=1 (all trial_heldout), "
          f"{advanced} at Phase>=3; both folds re-derive")


def check_leakage(j):
    """See REPRODUCIBILITY.md."""
    cmp_ = pd.read_csv(fd("leakage_controlled_comparison.csv"))
    # Loops both arms: the standard arm's ns marker is also on a published figure, needing the same check.
    expect_pos = {"standard_taskA": 727, "leakage_controlled": 535}
    for arm_name, a in j["arms"].items():
        assert a["n_pos"] == expect_pos[arm_name], (
            f"the {arm_name} arm is {expect_pos[arm_name]} positives, JSON says {a['n_pos']}")
        rows = cmp_[cmp_["evaluation"] == arm_name].set_index("method")["auc"]
        for stat, method in [("auc_pu", "PU full model"), ("auc_gps", "GPS overall")]:
            assert round(a[stat], 3) == rows[method], (
                f"leakage_controlled_delong.json {arm_name}.{stat}={a[stat]:.4f} rounds to "
                f"{round(a[stat], 3)}, but the comparison sheet says {rows[method]} for {method!r}")
            assert a["n_pos"] == int(
                cmp_.loc[(cmp_["method"] == method) & (cmp_["evaluation"] == arm_name),
                         "n_pos"].iloc[0]), \
                f"n_pos disagrees between the DeLong JSON and the {method!r} {arm_name} row"
        assert abs(a["delta"] - (a["auc_pu"] - a["auc_gps"])) < 1e-15, \
            f"{arm_name} delta {a['delta']} is not auc_pu - auc_gps"
        # See REPRODUCIBILITY.md.
        p = a["p_holm"]
        want = ("****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2
                else "*" if p < 5e-2 else "ns")
        assert a["stars"] == want, \
            f"{arm_name} stars={a['stars']!r} does not match p_holm={p:.4g} under the repo thresholds"
        print(f"  leakage DeLong [{arm_name}]: {a['auc_pu']:.4f} vs {a['auc_gps']:.4f} on "
              f"{a['n_pos']} positives, p = {p:.3g} [{a['stars']}]; both round to the sheet")
    assert len(cmp_) == 6, f"the comparison sheet is 6 rows since L2G was dropped, got {len(cmp_)}"
    assert "L2G" not in set(cmp_["method"]), "L2G was removed from this analysis; see the producer"


def check_vignette_knn(j):
    """See REPRODUCIBILITY.md."""
    assert j["condition"] == "Stim48hr", \
        f"the published neighborhood is Stim48hr only, this JSON says {j['condition']!r}"
    assert j["focal"] == "STAT4", f"the vignette gene is STAT4, this JSON says {j['focal']!r}"
    assert j["n_anchor_rows"] == 112, \
        f"expected 112 approved targets with a signature, got {j['n_anchor_rows']}"
    assert j["n_unlabeled_rows"] == 7044, \
        f"expected 7,044 unlabeled genes with a signature, got {j['n_unlabeled_rows']}"

    # build_vignette_meta.py copies these into vignette_meta.json, read by the renderer; never checked against source.
    meta = load_json("vignette_meta.json")
    for k in ("knn_p95", "knn_median_unlabeled"):
        assert meta[k] == j[k], \
            f"{k} is {j[k]} in vignette_knn_stats.json and {meta[k]} in vignette_meta.json"
    print(f"  knn: {j['n_anchor_rows']} anchors vs {j['n_unlabeled_rows']:,} unlabeled at "
          f"{j['condition']}; p95 and median match vignette_meta.json")


def main():
    heldout = load_json("heldout_trial_auc.json")
    groups = load_json("score_by_group_stats.json")
    scrambled = load_json("scrambled_feature_control.json")
    stages = load_json("immune_stages_counts.json")
    leakage = load_json("leakage_controlled_delong.json")
    knn = load_json("vignette_knn_stats.json")

    check_heldout(heldout)
    check_score_by_group(groups)
    check_scrambled(scrambled)
    check_permutation(scrambled)
    check_immune_stages(stages)
    check_leakage(leakage)
    check_vignette_knn(knn)

    out = {
        "heldout_trial_stats.csv": heldout_trial_stats(heldout),
        "score_by_group_summary.csv": score_by_group_summary(groups),
        "score_by_group_tests.csv": score_by_group_tests(groups),
        "scrambled_control_stats.csv": scrambled_control_stats(scrambled),
        "immune_stages_bootstrap.csv": immune_stages_bootstrap(stages),
        "immune_stages_depth_sweep.csv": immune_stages_depth_sweep(stages),
        "leakage_controlled_delong_stats.csv": leakage_controlled_delong_stats(leakage),
        "vignette_knn_summary.csv": vignette_knn_summary(knn),
    }
    for name, df in out.items():
        df.to_csv(fd(name), index=False)
        print(f"wrote figure_data/{name}  ({len(df)} rows x {len(df.columns)} cols)")


if __name__ == "__main__":
    main()
