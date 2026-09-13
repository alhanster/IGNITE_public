#!/usr/bin/env python
"""Backing tables for the genetics-vs-full-model discordance figure (S2).

Both arms score every gene; scores are converted to rank percentiles over the same
universe and combined as

    rank_shift = pctile_full - pctile_genetics

Positive values mark genes promoted by functional genomics, negative values
mark genes the genetic prior alone ranked higher.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import (chi2, false_discovery_control, fisher_exact,
                         mannwhitneyu, spearmanr)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
FIGDATA = os.path.join(ROOT, "figure_data")
SCRATCH = os.path.join(ROOT, "outputs", "discordance")

# Importing pu_target_model enforces the xgboost pin via its version guard; this file only calls pum.pu_bag.
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "model"))
import pu_target_model as pum  # noqa: E402
import _platform  # noqa: E402

GENETIC = ["lof.oe_ci.upper", "mis.z_score", "IEI", "gwas_score", "gene_burden_score"]
TAIL = 0.05                       # tail fraction for the druggability test
HI, LO = 0.75, 0.25               # core-set cuts: promoted = pctile_genetics < LO and pctile_full >= HI; demoted is the mirror
STRATA = [0, 0.50, 0.80, 0.95, 1.0]
STRATA_LAB = ["<50", "50-80", "80-95", ">95"]
N_NULL = 20                       # matched null sets for the FEATURE test only; the
                                  # Network null distribution size is 500, set in the fetch script; reported under stats key n_null_network.

DRIVER_FEATURES = [
    ("expected_n_regulators_residuals", "Regulatory burden\n(residual)"),
    ("reg_burden_sig_Stim48hr", "Regulator burden\n(48 h stim)"),
    ("zscore_Th2", "Th2 DE z-score"),
    ("zscore_Th17", "Th17 DE z-score"),
    ("has_cytokine", "Cytokine measured"),
    ("mis.z_score", "Missense constraint"),
    ("lof.oe_ci.upper", "LoF tolerance"),
    ("gwas_score", "GWAS score"),
    ("IEI", "IEI gene"),
    ("gene_burden_score", "Rare-variant burden"),
]

# See PROVENANCE.md.
CACHE_NET = os.path.join(FIGDATA, "discordance_network.csv")
CACHE_NODES = os.path.join(FIGDATA, "discordance_network_nodes.csv")
CACHE_GO = os.path.join(FIGDATA, "discordance_go_enrichment.csv")
CACHE_DARK = os.path.join(FIGDATA, "discordance_dark_proteome.csv")

REFRESH_HINT = (
    "  Refresh the cache, in this order:\n"
    "    1. make tables               (aborts here, having written the new core sets)\n"
    "    2. make discordance-network  (opt-in; ~10 min, one-time 135 MB download)\n"
    "    3. make tables               (passes)")


def _json_num(x):
    """Coerces non-finite values to None: NaN/Infinity from discordance_dark_proteome.csv would otherwise break jsonlite::fromJSON for the whole file."""
    v = float(x)
    return None if not np.isfinite(v) else v


def pctile(x):
    """Rank percentile over the full gene universe; ties are averaged."""
    return pd.Series(x).rank(pct=True).to_numpy()


def mantel_haenszel(tables):
    """Pooled 2x2 association across strata; each table entry is an (a, b, c, d) cell tuple in that order."""
    num = den = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        num += a - (a + b) * (a + c) / n
        den += (a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1))
    stat = (abs(num) - 0.5) ** 2 / den if den > 0 else np.nan
    return stat, float(chi2.sf(stat, 1))


# Both arms are fit once, at the production configuration.
m = pd.read_parquet(pum.MATRIX)
FEATS = [c for c in m.columns if c not in pum.META]
FG = [c for c in FEATS if c not in GENETIC]
role = m["pu_role"].to_numpy()
P_idx = np.where(role == "P")[0]
U_idx = np.where(role == "unlabeled")[0]
Xg = m[GENETIC].apply(pd.to_numeric, errors="coerce").to_numpy(float)
Xf = m[FEATS].apply(pd.to_numeric, errors="coerce").to_numpy(float)

pg = pctile(pum.pu_bag(Xg, P_idx, U_idx, GENETIC, T=pum.T_BAG, seed=pum.SEED)[0])
pf = pctile(pum.pu_bag(Xf, P_idx, U_idx, FEATS, T=pum.T_BAG, seed=pum.SEED)[0])
print("  fitted both arms at seed %d, T=%d" % (pum.SEED, pum.T_BAG), flush=True)

# Core sets from the quadrant rule are mutually exclusive by construction (LO < HI); no tie-break is needed.
is_up = (pg < LO) & (pf >= HI)
is_dn = (pg >= HI) & (pf < LO)

d = m[list(pum.META)].copy()
d["pctile_genetics"] = pg
d["pctile_full"] = pf
d["rank_shift"] = pf - pg
d["core"] = np.select([is_up, is_dn], ["promoted", "demoted"], default="none")
d["is_approved"] = (d.furthest_stage == "Approved").astype(int)
d["is_trial_heldout"] = (d.pu_role == "trial_heldout").astype(int)
for c in FEATS:
    d[c] = pd.to_numeric(m[c], errors="coerce")

# Annotations that never entered the training set.
drug = set(pd.read_csv(os.path.join(ROOT, "data", "reference",
                                    "druggable_genome_gene_list.csv"))
           .gene.str.strip())
d["is_druggable"] = d.gene.isin(drug).astype(int)

# Matching covariates are computed here so they land as columns 6-7 of discordance_scatter.csv; this consumes no RNG draws.
obs_cols = ["zscore_Th1", "zscore_Th2", "zscore_Th17", "zscore_Treg"]
cau_cols = [c for c in FEATS if c.startswith(("expected_n", "polar_", "reg_burden",
                                              "n_sig_"))]
d["cov_class"] = np.select(
    [d[obs_cols].notna().sum(axis=1).gt(0) & d[cau_cols].notna().sum(axis=1).gt(0),
     d[obs_cols].notna().sum(axis=1).gt(0),
     d[cau_cols].notna().sum(axis=1).gt(0)],
    ["both", "obs", "cau"], default="none")
d["gen_dec"] = pd.qcut(d.pctile_genetics, 10, labels=False)

os.makedirs(FIGDATA, exist_ok=True)
os.makedirs(SCRATCH, exist_ok=True)

# Column order matters: the renderer reads the first five columns and the fetch target reads the last two.
d[["gene", "core", "pctile_genetics", "pctile_full", "rank_shift",
   "cov_class", "gen_dec"]].to_csv(
    os.path.join(FIGDATA, "discordance_scatter.csv"), index=False)

# Drivers: no threshold.
rows = []
for col, label in DRIVER_FEATURES:
    rho, p = spearmanr(d[col], d.rank_shift, nan_policy="omit")
    rows.append((col, label, rho, p, int(d[col].notna().sum())))
pd.DataFrame(rows, columns=["feature", "label", "rho", "p", "n_nonnull"]) \
  .sort_values("rho").to_csv(os.path.join(FIGDATA, "discordance_drivers.csv"),
                             index=False)

# Druggability tails: exploratory; these tails are not the marked gene set.
hi = d[d.rank_shift >= d.rank_shift.quantile(1 - TAIL)]
lo = d[d.rank_shift <= d.rank_shift.quantile(TAIL)]
rows = []
for lab, sel in [("promoted", hi), ("demoted", lo)]:
    rest = d.drop(sel.index)
    k, n = int(sel.is_druggable.sum()), len(sel)
    K, N = int(rest.is_druggable.sum()), len(rest)
    odds, p = fisher_exact([[k, n - k], [K, N - K]])
    rows.append((lab, n, k, 100 * k / n, 100 * K / N, odds, p))
pd.DataFrame(rows, columns=["direction", "n", "k_druggable", "pct_sel",
                            "pct_rest", "odds_ratio", "p"]) \
  .to_csv(os.path.join(SCRATCH, "druggability.csv"), index=False)

# Validation: stratified median split.
d["gen_stratum"] = pd.cut(d.pctile_genetics, STRATA, labels=STRATA_LAB)
rows, tables = [], []
for lab in STRATA_LAB:
    s = d[d.gen_stratum == lab]
    med = s.rank_shift.median()
    up, dn = s[s.rank_shift > med], s[s.rank_shift <= med]
    a, c = int(up.is_trial_heldout.sum()), int(dn.is_trial_heldout.sum())
    b, e = len(up) - a, len(dn) - c
    odds, p = fisher_exact([[a, b], [c, e]])
    rows.append((lab, len(s), 100 * a / len(up), 100 * c / len(dn), a, c, odds, p))
    tables.append((a, b, c, e))
val = pd.DataFrame(rows, columns=["gen_stratum", "n", "pct_up", "pct_down",
                                  "k_up", "k_down", "odds_ratio", "fisher_p"])
val.to_csv(os.path.join(FIGDATA, "discordance_validation.csv"), index=False)
mh_stat, mh_p = mantel_haenszel(tables)

# Evidence independence: tests whether IEI/GWAS evidence already captures the core sets defined by the quadrant rule.
groups = {
    "core promoted": d[d.core == "promoted"],
    "core demoted": d[d.core == "demoted"],
    "concordant high": d[(d.pctile_genetics >= HI) & (d.pctile_full >= HI)],
    "all genes": d,
}
rows = []
for lab, s in groups.items():
    rows.append((lab, len(s), int(s.IEI.fillna(0).sum()),
                 100 * s.IEI.fillna(0).mean(),
                 100 * (s.gwas_score.fillna(0) > 0).mean(),
                 100 * s.is_druggable.mean(),
                 float(s["mis.z_score"].median())))
pd.DataFrame(rows, columns=["group", "n", "k_IEI", "pct_IEI", "pct_gwas_nonzero",
                            "pct_druggable", "median_mis_z"]) \
  .to_csv(os.path.join(FIGDATA, "discordance_evidence.csv"), index=False)

# Feature drivers: marked genes vs controls matched on FG-coverage class x genetics decile, pooling N_NULL draws into one control set.
rng = np.random.default_rng(pum.SEED)
pool = d[d.core == "none"]


def matched_controls(genes, n_sets=N_NULL):
    """n_sets random gene sets matching `genes` jointly on coverage class and genetics-score decile."""
    tgt, out = d[d.gene.isin(genes)], []
    for _ in range(n_sets):
        picks = []
        for (cc, gd), grp in tgt.groupby(["cov_class", "gen_dec"], observed=True):
            cand = pool[(pool.cov_class == cc) & (pool.gen_dec == gd)]
            if len(cand):
                picks += list(cand.sample(n=min(len(grp), len(cand)),
                                          random_state=int(rng.integers(1e6))).gene)
        out.append(sorted(picks))
    return out


rows = []
for lab in ("promoted", "demoted"):
    genes = sorted(d.loc[d.core == lab, "gene"])
    ctrl_sets = matched_controls(genes)
    ctrl = d[d.gene.isin({g for s_ in ctrl_sets for g in s_})]
    sel = d[d.gene.isin(genes)]
    for f in FEATS:
        a = pd.to_numeric(sel[f], errors="coerce").dropna()
        b = pd.to_numeric(ctrl[f], errors="coerce").dropna()
        # A constant feature has no defined U statistic; IEI is constant (all 0) among promoted genes, itself a finding (see evidence.csv).
        if len(a) < 5 or len(b) < 20 or a.nunique() <= 1:
            continue
        u, p = mannwhitneyu(a, b)
        rows.append((lab, f, "genetic" if f in GENETIC else "functional genomics",
                     len(a), float(a.median()), float(b.median()),
                     2 * u / (len(a) * len(b)) - 1, p))
feat = pd.DataFrame(rows, columns=["direction", "feature", "block", "n_sel",
                                   "med_sel", "med_ctrl", "rank_biserial", "p"])
feat["fdr"] = np.nan
for lab in feat.direction.unique():
    msk = feat.direction == lab
    feat.loc[msk, "fdr"] = false_discovery_control(
        feat.loc[msk, "p"].to_numpy(float), method="bh")
feat.sort_values(["direction", "p"]).to_csv(
    os.path.join(FIGDATA, "discordance_features.csv"), index=False)

core = d[d.core != "none"][
    ["gene", "core", "pctile_genetics",
     "pctile_full", "rank_shift", "IEI", "gwas_score", "mis.z_score",
     "lof.oe_ci.upper", "expected_n_regulators_residuals", "is_druggable",
     "is_approved", "is_trial_heldout"]
].sort_values(["core", "rank_shift"])
core.to_csv(os.path.join(FIGDATA, "discordance_core_genes.csv"), index=False)

# See REPRODUCIBILITY.md.
expect = {"promoted": int((d.core == "promoted").sum()),
          "demoted": int((d.core == "demoted").sum())}
core_genes = set(d.loc[d.core != "none", "gene"])
net_stats = {}
stats_dark = {}

# See REPRODUCIBILITY.md.
problems = []
n_dark_prom = None

if os.path.exists(CACHE_NET):
    net = pd.read_csv(CACHE_NET)
    got = {r.direction: int(r.n_genes) for r in net.itertuples()}
    for direction, n in expect.items():
        if got.get(direction) != n:
            problems.append("network cache says %s=%s, current core set is %d"
                            % (direction, got.get(direction), n))

if os.path.exists(CACHE_NODES):
    nod = pd.read_csv(CACHE_NODES)
    orphan = sorted(set(nod.gene) - core_genes)
    if orphan:
        problems.append("%d cached network node(s) are no longer core genes: %s"
                        % (len(orphan), orphan[:8]))

if os.path.exists(CACHE_DARK):
    # Uses .iloc[0] on a filtered frame rather than .loc, so a duplicated group label raises this block's own error rather than a generic pandas one.
    dk = pd.read_csv(CACHE_DARK)
    prom_rows = dk[dk.group == "core promoted"] if "group" in dk.columns else dk.iloc[:0]
    if len(prom_rows) and "n" in dk.columns and "n_dark" in dk.columns:
        if int(prom_rows.iloc[0]["n"]) != expect["promoted"]:
            problems.append("dark-proteome cache says promoted n=%d, current is %d"
                            % (int(prom_rows.iloc[0]["n"]), expect["promoted"]))
        n_dark_prom = int(prom_rows.iloc[0]["n_dark"])

# See REPRODUCIBILITY.md.
if os.path.exists(CACHE_GO):
    go_ = pd.read_csv(CACHE_GO)
    gprom = go_[go_.direction == "promoted"] if "direction" in go_.columns else go_.iloc[:0]
    if len(gprom) and "n_fg" in go_.columns:
        if n_dark_prom is None:
            problems.append(
                "GO cache is present but the dark-proteome cache is missing or unreadable, "
                "so n_fg cannot be checked against the core set -- the two are only "
                "meaningful together")
        else:
            n_fg = int(gprom.iloc[0]["n_fg"])
            if n_fg + n_dark_prom != expect["promoted"]:
                problems.append(
                    "GO cache says n_fg=%d annotated promoted genes and the dark-proteome "
                    "cache says %d unannotated, totalling %d -- but the current promoted "
                    "set is %d" % (n_fg, n_dark_prom, n_fg + n_dark_prom,
                                   expect["promoted"]))

if problems and not _platform.on_reference_platform():
    # The cache describes the reference platform's core sets; a refit elsewhere can move them by a
    # gene. Reported, not fatal, off the reference platform. See REPRODUCIBILITY.md.
    print("  WARNING (not the reference platform): STRING/GO cache describes different core sets:\n"
          + "".join("    %s\n" % p for p in problems))
    problems = []
if problems:
    raise SystemExit(
        "\nSTRING/GO CACHE IS STALE -- the core sets moved since it was built.\n"
        + "".join("    %s\n" % p for p in problems)
        + "\n  These tables are committed inputs, not rebuilt by `make tables`:\n"
          "    figure_data/discordance_{network,network_edges,network_nodes}.csv\n"
          "    figure_data/discordance_go_*.csv\n"
          "    figure_data/discordance_dark_proteome*.csv\n\n"
        + REFRESH_HINT + "\n\n"
          "  discordance_stats.json was NOT written, so nothing downstream can\n"
          "  quote a p-value computed on a different gene set.\n")

if os.path.exists(CACHE_NET):
    for r in net.itertuples():
        net_stats["net_%s_edges" % r.direction] = int(r.n_edges)
        net_stats["net_%s_null_mean" % r.direction] = float(r.null_mean_edges)
        net_stats["net_%s_p" % r.direction] = float(r.p_empirical)
        net_stats["net_%s_connected" % r.direction] = int(r.n_connected)
    # The network null's draw count is a separate key from N_NULL, so a missing cache cannot silently mislabel panel b's caption with the feature test's draw count.
    if "n_null_draws" in net.columns:
        net_stats["n_null_network"] = int(net.n_null_draws.iloc[0])
else:
    print("   note: no %s -- panel b's stats keys are omitted; run "
          "`make discordance-network`" % os.path.relpath(CACHE_NET, ROOT))

if os.path.exists(CACHE_GO):
    g_ = pd.read_csv(CACHE_GO)
    gp = g_[(g_.direction == "promoted") & (g_.fdr < 0.05)].sort_values("p")
    net_stats["go_n_sig_promoted"] = int(len(gp))
    net_stats["go_n_sig_demoted"] = int(
        ((g_.direction == "demoted") & (g_.fdr < 0.05)).sum())
    # Keyed on the top promoted term by p, not on the significant slice. The renderer stopped
    # gating these keys when S2's GO panel was dropped, so a run with nothing under FDR 0.05 would
    # otherwise omit four keys silently instead of reporting a top term with a weak FDR.
    top = g_[g_.direction == "promoted"].sort_values("p")
    if len(top):
        net_stats["go_top_term"] = str(top.iloc[0]["term"])
        net_stats["go_top_or"] = float(top.iloc[0]["odds_ratio"])
        net_stats["go_top_fdr"] = float(top.iloc[0]["fdr"])
        net_stats["go_n_background"] = int(top.iloc[0]["n_bg"])
else:
    print("   note: no %s -- the GO stats keys are omitted; run "
          "`make discordance-network`" % os.path.relpath(CACHE_GO, ROOT))

# Dark proteome (from the fetch target's GO-annotation GAF): dark means no GO biological-process term, not uncharacterised protein; UniProt annotation depth shows no such difference.
if os.path.exists(CACHE_DARK):
    dk = pd.read_csv(CACHE_DARK)
    prom = dk[dk.group == "core promoted"]
    cm = dk[dk.group.str.startswith("coverage-matched")]
    if len(prom) and len(cm):
        stats_dark = {
            "dark_promoted": int(prom.iloc[0]["n_dark"]),
            "dark_promoted_pct": _json_num(prom.iloc[0]["pct_dark"]),
            "dark_background_pct": _json_num(cm.iloc[0]["pct_dark"]),
            "dark_OR": _json_num(cm.iloc[0]["OR_vs_promoted"]),
            "dark_p": _json_num(cm.iloc[0]["fisher_p"]),
        }
else:
    print("   note: no %s -- the dark-proteome keys are omitted; run "
          "`make discordance-network`" % os.path.relpath(CACHE_DARK, ROOT))

# Provenance check: rho < 1 means the refit arm's seed, T_BAG, tree depth or xgboost build does not match the committed model, so the run aborts.
committed = os.path.join(pum.EVAL, "full_model_pu_scores.csv")
if not os.path.exists(committed):
    raise SystemExit(
        "\nMissing %s\n  Run `make tables` from the start -- pu_target_model.py "
        "must run before this step.\n" % os.path.relpath(committed, ROOT))
g = pd.read_csv(committed)
col = next((c for c in ("score_full", "pu_score", "score") if c in g.columns), None)
if col is None:
    raise SystemExit("\n%s has no recognisable score column (looked for "
                     "score_full / pu_score / score)\n"
                     % os.path.relpath(committed, ROOT))
j = d[["gene", "pctile_full"]].merge(g[["gene", col]], on="gene")
rho_committed = float(spearmanr(j.pctile_full, j[col])[0])
if not np.isclose(rho_committed, 1.0, atol=1e-9):
    raise SystemExit(
        "\nFULL ARM DOES NOT REPRODUCE THE COMMITTED RANKING\n"
        "    spearman_vs_committed_full = %.10f, expected 1.0 (n=%d genes joined)\n\n"
        "  The arm refit here is not the model behind outputs/model/full_model_pu_scores.csv,\n"
        "  so every gene this analysis names is suspect. Usual causes: a changed\n"
        "  pum.SEED / pum.T_BAG / tree depth, or an xgboost other than the pin.\n"
        "  Check `make check-versions` first.\n" % (rho_committed, len(j)))

stats = {
    "seed": int(pum.SEED),
    "bags": int(pum.T_BAG),
    # Tree depth is read from the fitted learner, not the MODEL_MAX_DEPTH environment variable, so the recorded depth matches the depth actually used.
    "max_depth": int(pum.base_learner(0).max_depth),
    "n_genes": int(len(d)),
    "n_features_full": len(FEATS),
    "n_features_fg": len(FG),
    "hi_cut": HI,
    "lo_cut": LO,
    "tail_fraction": TAIL,
    "spearman_arms": float(spearmanr(d.pctile_genetics, d.pctile_full)[0]),
    "n_core_promoted": expect["promoted"],
    "n_core_demoted": expect["demoted"],
    "mh_chi2": float(mh_stat),
    "mh_p": float(mh_p),
    "n_strata_up_gt_down": int((val.pct_up > val.pct_down).sum()),
    "n_strata": len(val),
    "spearman_vs_committed_full": rho_committed,
    "n_null_features": N_NULL,
}
stats.update(net_stats)
stats.update(stats_dark)
with open(os.path.join(FIGDATA, "discordance_stats.json"), "w") as fh:
    json.dump(stats, fh, indent=2)

print("wrote", os.path.relpath(os.path.join(FIGDATA, "discordance_stats.json"), ROOT))
print("   seed %d, T=%d, depth=%d | arms rho %.4f"
      % (pum.SEED, pum.T_BAG, stats["max_depth"], stats["spearman_arms"]))
print("   core: %d promoted / %d demoted (quadrant rule at %.2f/%.2f)"
      % (expect["promoted"], expect["demoted"], HI, LO))
print("   MH chi2=%.2f p=%.2g (%d/%d strata) | rho vs committed full = %.4f"
      % (mh_stat, mh_p, stats["n_strata_up_gt_down"], stats["n_strata"],
         rho_committed))
