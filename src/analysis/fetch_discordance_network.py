#!/usr/bin/env python
"""STRING coherence, GO over-representation and dark-proteome tests for the discordance figure's core sets. Opt-in, run via make discordance-network.

See REPRODUCIBILITY.md, discordance-network, for archive provenance, rerun conditions and the null-model rationale.
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import urllib.request
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
assert os.path.isfile(os.path.join(ROOT, "Makefile")), f"repo root not found: {ROOT}"
FIGDATA = os.path.join(ROOT, "figure_data")
RAW = os.path.join(ROOT, "data", "raw", "discordance")

MIN_SCORE = 400          # STRING medium confidence (0-1000 scale in the flat file)
SEED = 4                 # cross-checked against discordance_stats.json below

STRING_BASE = "https://stringdb-downloads.org/download"
STRING_FILES = {
    "links": "%s/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" % STRING_BASE,
    "info": "%s/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" % STRING_BASE,
    "aliases": ("%s/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
                % STRING_BASE),
}
GAF_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz"
UNIPROT_URL = ("https://rest.uniprot.org/uniprotkb/stream"
               "?query=organism_id:9606+AND+reviewed:true"
               "&fields=gene_primary,annotation_score,cc_function"
               "&format=tsv&compressed=true")
QUICKGO = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%s"

EDGE_CACHE = os.path.join(RAW, "string_edges_ge%d.csv" % MIN_SCORE)
# Sidecar recording the gene universe the edge cache was derived over. See REPRODUCIBILITY.md, What a green gate does not cover.
EDGE_META = os.path.join(RAW, "string_edges_ge%d.meta.json" % MIN_SCORE)
GAF_CACHE = os.path.join(RAW, "goa_human.gaf.gz")
UNIPROT_CACHE = os.path.join(RAW, "uniprot_human_reviewed.tsv.gz")
LABEL_CACHE = os.path.join(ROOT, "data", "reference", "go_term_labels.json")

# GOA GAF and UniProt endpoints are unpinned; STRING is pinned to v12.0. See REPRODUCIBILITY.md, discordance-network.


def cached(url, name):
    """Downloads into data/raw/discordance/ once; reused on subsequent runs."""
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, name)
    if not os.path.exists(path):
        print("downloading %s ..." % name, flush=True)
        urllib.request.urlretrieve(url, path)
    return path


def need(name):
    p = os.path.join(FIGDATA, name)
    if not os.path.exists(p):
        raise SystemExit("\nMissing %s\n  Run `make tables` first.\n"
                         % os.path.relpath(p, ROOT))
    return p


# ---- committed core sets ----
d = pd.read_csv(need("discordance_scatter.csv"))
S = json.load(open(need("discordance_stats.json")))
HI = float(S["hi_cut"])
if int(S["seed"]) != SEED:
    raise SystemExit(
        "\nSEED MISMATCH: this script uses %d, discordance_stats.json records %d.\n"
        "  The matched nulls below would be drawn under a different seed than the\n"
        "  core sets they are matched to. Reconcile SEED here with pum.SEED.\n"
        % (SEED, int(S["seed"])))
for col in ("cov_class", "gen_dec"):
    if col not in d.columns:
        raise SystemExit(
            "\ndiscordance_scatter.csv has no `%s` column.\n"
            "  It is written by build_discordance.py for exactly this script; a\n"
            "  scatter table without it predates that and must be rebuilt:\n"
            "    make tables\n" % col)


def parse_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=500,
                    help="matched null draws per direction (default 500)")
    ap.add_argument("--validate", action="store_true",
                    help="check the offline edge counter against the committed "
                         "cache, then exit")
    return ap.parse_args(argv)


args = parse_args(sys.argv[1:])
N_NULL = args.n_null

rng = np.random.default_rng(SEED)
pool = d[d.core == "none"]


def matched(genes, n_sets=N_NULL):
    """Returns n_sets gene sets matching `genes` on coverage class x genetics decile."""
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


# Offline STRING edge table, built once and cached. See PROVENANCE.md, Written by the opt-in `make discordance-network` target (6 files).

def _symbol_to_ensp(universe):
    """Maps preferred_name to ENSP for symbols in `universe`, using aliases only when the preferred name does not resolve."""
    m = {}
    with gzip.open(cached(STRING_FILES["info"],
                          "9606.protein.info.v12.0.txt.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1] in universe:
                m.setdefault(f[1], f[0])
    unmapped = universe - set(m)
    if unmapped:
        # Loads the 20 MB alias file only as a fallback when preferred_name lookup misses.
        print("   %d symbol(s) not a STRING preferred_name; trying aliases"
              % len(unmapped), flush=True)
        with gzip.open(cached(STRING_FILES["aliases"],
                              "9606.protein.aliases.v12.0.txt.gz"), "rt") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1] in unmapped and f[1] not in m:
                    m[f[1]] = f[0]
    return m


def _universe_fingerprint(universe):
    """Produces a stable id for the gene set an edge cache was derived over."""
    genes = sorted(set(map(str, universe)))
    return {"min_score": MIN_SCORE, "n_universe": len(genes),
            "universe_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest()}


def load_edges(universe):
    """Returns edges with score >= MIN_SCORE, restricted to gene pairs within `universe`."""
    want = _universe_fingerprint(universe)
    if os.path.exists(EDGE_CACHE):
        # Aborts on a cache miss rather than re-deriving from the 135 MB STRING download, to avoid a confusing offline failure or a quietly wrong number.
        got = json.load(open(EDGE_META)) if os.path.exists(EDGE_META) else None
        if got != want:
            raise SystemExit(
                "STRING edge cache was built over a different gene universe.\n"
                "  cache:    %s\n"
                "  expected: %s\n"
                "  found:    %s\n"
                "Every edge touching a gene absent from the cached universe is missing,\n"
                "so the network p-value would be computed on a truncated graph.\n"
                "  rm %s %s\n"
                "then rerun `make discordance-network` to re-derive it."
                % (os.path.relpath(EDGE_CACHE, ROOT),
                   want, got if got is not None else "no sidecar (cache predates this check)",
                   os.path.relpath(EDGE_CACHE, ROOT), os.path.relpath(EDGE_META, ROOT)))
        e = pd.read_csv(EDGE_CACHE)
        return {frozenset((a, b)): s for a, b, s in zip(e.a, e.b, e.score)}
    sym = _symbol_to_ensp(set(universe))
    ensp2sym = {v: k for k, v in sym.items()}
    edges = {}
    with gzip.open(cached(STRING_FILES["links"],
                          "9606.protein.links.v12.0.txt.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            p1, p2, sc = line.split()
            if int(sc) < MIN_SCORE:
                continue
            a, b = ensp2sym.get(p1), ensp2sym.get(p2)
            if a is None or b is None or a == b:
                continue
            # The STRING file lists each edge in both directions; keying by a frozenset of the pair dedupes them.
            edges[frozenset((a, b))] = int(sc) / 1000.0
    out = sorted((min(k), max(k), v) for k, v in edges.items())
    os.makedirs(RAW, exist_ok=True)
    pd.DataFrame(out, columns=["a", "b", "score"]).to_csv(EDGE_CACHE, index=False)
    with open(EDGE_META, "w") as fh:
        json.dump(want, fh, indent=1)
    print("   cached %d edges (>=%d) over %d genes -> %s"
          % (len(out), MIN_SCORE, want["n_universe"],
             os.path.relpath(EDGE_CACHE, ROOT)), flush=True)
    return edges


EDGES = load_edges(d.gene.dropna().unique())
MAPPED = {g for k in EDGES for g in k}


def network(symbols):
    """Returns the induced subgraph over `symbols` as (n_nodes, n_edges, n_connected, edges)."""
    s = set(symbols)
    hits = sorted(((min(k), max(k), v) for k, v in EDGES.items() if k <= s),
                  key=lambda r: (r[0], r[1]))
    deg = Counter()
    for a, b, _ in hits:
        deg[a] += 1
        deg[b] += 1
    edges = [{"a": a, "b": b, "score": v} for a, b, v in hits]
    return len(s), len(hits), len(deg), edges


def validate_against_cache():
    """Checks whether the offline counter reproduces the committed network cache. See REPRODUCIBILITY.md, What a green gate does not cover."""
    nod = pd.read_csv(need("discordance_network_nodes.csv"))
    ref = pd.read_csv(need("discordance_network_edges.csv"))
    net = pd.read_csv(need("discordance_network.csv")).set_index("direction")
    ok = True
    for direction in ("promoted", "demoted"):
        genes = sorted(nod.loc[nod.direction == direction, "gene"])
        _, n_e, n_c, edges = network(genes)
        exp_e = int(net.loc[direction, "n_edges"])
        exp_c = int(net.loc[direction, "n_connected"])
        r = ref[ref.direction == direction]
        got = sorted((min(e["a"], e["b"]), max(e["a"], e["b"]), round(e["score"], 3))
                     for e in edges)
        want = sorted((min(a, b), max(a, b), round(s, 3))
                      for a, b, s in zip(r.a, r.b, r.score))
        match = (n_e, n_c, got) == (exp_e, exp_c, want)
        ok &= match
        print("%-9s n=%2d | recomputed %2d edges/%2d connected | committed %2d/%2d | %s"
              % (direction, len(genes), n_e, n_c, exp_e, exp_c,
                 "MATCH" if match else "MISMATCH"))
        if got != want:
            for lab, s in (("recomputed only", set(got) - set(want)),
                           ("committed only", set(want) - set(got))):
                if s:
                    print("     %s: %s" % (lab, sorted(s)))
    print("\nvalidation %s" % (
        "PASSED -- the offline counter reproduces the committed cache" if ok
        else "FAILED -- the STRING source has moved; do NOT overwrite the cache"))
    return ok


if args.validate:
    sys.exit(0 if validate_against_cache() else 1)

os.makedirs(FIGDATA, exist_ok=True)

rows, all_edges = [], []
for direction in ("promoted", "demoted"):
    genes = sorted(d.loc[d.core == direction, "gene"])
    n_nodes, n_edges, n_conn, edges = network(genes)
    n_mapped = sum(1 for g in genes if g in MAPPED)
    nulls = [network(s)[1] for s in matched(genes)]
    p = (sum(1 for x in nulls if x >= n_edges) + 1) / (len(nulls) + 1)
    rows.append((direction, n_nodes, n_mapped, n_edges, n_conn,
                 float(np.mean(nulls)), int(max(nulls)),
                 n_edges / max(np.mean(nulls), 0.01), p, N_NULL))
    print("%-9s %2d genes (%d in STRING) | %2d edges (%d connected) | "
          "null mean %.1f max %d | p=%.4f"
          % (direction, n_nodes, n_mapped, n_edges, n_conn, np.mean(nulls),
             max(nulls), p))
    for e in edges:
        all_edges.append((direction, e["a"], e["b"], e["score"]))

pd.DataFrame(rows, columns=["direction", "n_genes", "n_mapped", "n_edges",
                            "n_connected", "null_mean_edges", "null_max_edges",
                            "fold", "p_empirical", "n_null_draws"]).to_csv(
    os.path.join(FIGDATA, "discordance_network.csv"), index=False)
ed = pd.DataFrame(all_edges, columns=["direction", "a", "b", "score"])
ed.to_csv(os.path.join(FIGDATA, "discordance_network_edges.csv"), index=False)

deg = {}
for _, r in ed.iterrows():
    for g in (r.a, r.b):
        deg[(r.direction, g)] = deg.get((r.direction, g), 0) + 1
pd.DataFrame(
    [(dr, g, deg.get((dr, g), 0))
     for dr in ("promoted", "demoted")
     for g in sorted(d.loc[d.core == dr, "gene"])],
    columns=["direction", "gene", "degree"]
).sort_values(["direction", "degree"], ascending=[True, False]).to_csv(
    os.path.join(FIGDATA, "discordance_network_nodes.csv"), index=False)

# ---- GO over-representation against the full non-discordant gene set as background ----
if not os.path.exists(GAF_CACHE):
    cached(GAF_URL, os.path.basename(GAF_CACHE))
if not os.path.exists(UNIPROT_CACHE):
    print("downloading UniProt annotation-depth table ...", flush=True)
    cached(UNIPROT_URL, os.path.basename(UNIPROT_CACHE))

pairs = set()
with gzip.open(GAF_CACHE, "rt") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 15 or f[8] != "P" or f[3].startswith("NOT"):
            continue
        pairs.add((f[2], f[4]))

g2t = {}
for sym, term in pairs:
    g2t.setdefault(sym, set()).add(term)


def ora(fg_genes, bg_genes, min_size=10, max_size=500, min_k_display=3):
    """Over-representation with a pre-specified correction family.

    The BH family is every GO term whose annotated universe size (foreground plus
    background) falls in [min_size, max_size] -- a filter on term size only, which is
    independent of the Fisher statistic under the null (Bourgon, Gentleman & Huber,
    PNAS 2010). min_k_display controls what is shown, never what is corrected.
    """
    fg = [g for g in fg_genes if g in g2t]
    bg = [g for g in bg_genes if g in g2t]
    cf, cb = Counter(), Counter()
    for g in fg:
        cf.update(g2t[g])
    for g in bg:
        cb.update(g2t[g])
    universe = {t: cf.get(t, 0) + cb.get(t, 0) for t in set(cf) | set(cb)}
    family = sorted(t for t, n in universe.items() if min_size <= n <= max_size)
    rows = []
    for term in family:
        k = cf.get(term, 0)
        K = cb.get(term, 0)
        odds, p = fisher_exact([[k, len(fg) - k], [K, len(bg) - K]],
                               alternative="greater")
        rows.append((term, k, len(fg), K, len(bg), odds, p))
    r = pd.DataFrame(rows, columns=["go_id", "k_fg", "n_fg", "k_bg", "n_bg",
                                    "odds_ratio", "p"])
    if len(r):
        r["fdr"] = false_discovery_control(r.p.to_numpy(float), method="bh")
        r["n_family"] = len(family)
        r["displayed"] = r.k_fg >= min_k_display
    else:
        # Without these an empty direction drops the columns, and the concat below turns
        # `displayed` to NaN for the other direction's rows.
        r = r.assign(fdr=pd.Series(dtype=float), n_family=pd.Series(dtype=int),
                     displayed=pd.Series(dtype=bool))
    return r.sort_values("p")


fg_sets = {k: sorted(d.loc[d.core == k, "gene"]) for k in ("promoted", "demoted")}
marked = set(fg_sets["promoted"]) | set(fg_sets["demoted"])
bg_full = sorted(set(d.gene) - marked)
cov_of = dict(zip(d.gene, d.cov_class))
bg_cov = sorted(g for g in bg_full if cov_of.get(g) == "both")

blocks = []
for direction, fg in fg_sets.items():
    for tag, bgset in (("full_universe", bg_full),
                       ("coverage_restricted", bg_cov)):
        r = ora(fg, bgset)
        if len(r):
            blocks.append(r.assign(direction=direction, background=tag))
go_all = pd.concat(blocks)

# Only resolved GO labels persist; caching an id as its own label leaves it unresolved forever.
labels = json.load(open(LABEL_CACHE)) if os.path.exists(LABEL_CACHE) else {}
labels = {k: v for k, v in labels.items() if v != k}      # drop past poisoning
missing = sorted(set(go_all.loc[go_all.displayed, "go_id"]) - set(labels))
resolved, unresolved = {}, []
for tid in missing:
    try:
        with urllib.request.urlopen(QUICKGO % tid, timeout=30) as r:
            name = json.load(r)["results"][0]["name"]
        resolved[tid] = name
        print("   resolved %s -> %s" % (tid, name))
    except Exception as exc:
        unresolved.append(tid)
        print("   WARNING: could not resolve %s (%s)" % (tid, exc))
if resolved:
    labels.update(resolved)
    with open(LABEL_CACHE, "w") as fh:
        json.dump(labels, fh, indent=2)
if unresolved:
    print("   WARNING: %d term(s) unresolved -- Supp Table 9's GO sheet will show raw GO ids for "
          "%s. Not cached, so the next run retries." % (len(unresolved), unresolved))
go_all["term"] = go_all.go_id.map(lambda t: labels.get(t, t))

cols = ["direction", "background", "go_id", "term", "k_fg", "n_fg", "k_bg",
        "n_bg", "odds_ratio", "p", "fdr", "n_family", "displayed"]
go_all[cols].to_csv(os.path.join(FIGDATA, "discordance_go_enrichment_full.csv"),
                    index=False)
# Already BH-corrected in ora(); re-correcting this slice alone would change every FDR in the primary table.
primary = go_all[(go_all.background == "coverage_restricted") & go_all.displayed].drop(
    columns=["background", "displayed"])
primary.to_csv(os.path.join(FIGDATA, "discordance_go_enrichment.csv"), index=False)

# Per-term diffs are shown explicitly since change direction is not uniform across terms.
sig_any = go_all[go_all.displayed][["direction", "term"]].drop_duplicates()
piv = (go_all.merge(sig_any, on=["direction", "term"])
       .pivot_table(index=["direction", "term"], columns="background",
                    values=["odds_ratio", "fdr"]).reset_index())
piv.columns = ["direction", "term", "fdr_cov", "fdr_full", "OR_cov", "OR_full"]
piv["pct_change_full_vs_cov"] = (100 * (piv.OR_full / piv.OR_cov - 1)).round(1)
piv["sig_flips"] = (piv.fdr_cov < 0.05) != (piv.fdr_full < 0.05)
piv.sort_values("pct_change_full_vs_cov", ascending=False).round(4).to_csv(
    os.path.join(FIGDATA, "discordance_go_background_comparison.csv"), index=False)
if piv.sig_flips.any():
    print("   note: %d term(s) cross FDR 0.05 between backgrounds: %s"
          % (int(piv.sig_flips.sum()), list(piv.loc[piv.sig_flips, "term"])))

print("GO: background %d coverage-matched non-discordant annotated genes"
      % len([g for g in bg_cov if g in g2t]))
for direction in ("promoted", "demoted"):
    n = int(((primary.direction == direction) & (primary.fdr < 0.05)).sum())
    print("  %-9s %d terms at FDR<0.05" % (direction, n))

# "Dark" means no GO BP term in the GAF; UniProt depth is unchanged, so it does not mean uncharacterised protein.
d["dark"] = (~d.gene.isin(set(g2t))).astype(int)

prom_mask = d.core == "promoted"
k_d, n_d = int(d.loc[prom_mask, "dark"].sum()), int(prom_mask.sum())
conc_mask = (d.pctile_genetics >= HI) & (d.pctile_full >= HI) & (d.core == "none")
groups = [
    ("core promoted", prom_mask),
    ("core demoted", d.core == "demoted"),
    ("concordant high (both arms >=%d%%)" % (HI * 100), conc_mask),
    ("coverage-matched non-discordant",
     (d.core == "none") & d.cov_class.isin(["both", "obs"])),
    ("all non-discordant genes", d.core == "none"),
]
dark_rows = []
for label, mask in groups:
    sub = d[mask]
    K, N = int(sub.dark.sum()), len(sub)
    if label == "core promoted" or N == 0:
        odds = pv = None
    else:
        odds, pv = fisher_exact([[k_d, n_d - k_d], [K, N - K]],
                                alternative="greater")
    # Checks `x is None`, not truthiness: 0.0 is a real odds ratio or p-value, not a missing comparison.
    dark_rows.append((label, N, K, None if N == 0 else round(100 * K / N, 1),
                      None if odds is None else round(odds, 3),
                      None if pv is None else round(pv, 6)))
pd.DataFrame(dark_rows, columns=["group", "n", "n_dark", "pct_dark",
                                 "OR_vs_promoted", "fisher_p"]).to_csv(
    os.path.join(FIGDATA, "discordance_dark_proteome.csv"), index=False)

# All three rows are computed live against the same coverage-matched background, so axes are comparable.
bg = d[(d.core == "none") & d.cov_class.isin(["both", "obs"])]
axes_rows = [(
    "GO biological process (no BP term)",
    "%.1f%% vs %.1f%%" % (100 * d.loc[prom_mask, "dark"].mean(),
                          100 * bg.dark.mean()),
    *[round(x, 4) for x in fisher_exact(
        [[k_d, n_d - k_d],
         [int(bg.dark.sum()), len(bg) - int(bg.dark.sum())]],
        alternative="greater")],
)]

u = pd.read_csv(UNIPROT_CACHE, sep="\t", low_memory=False)
u.columns = ["gene", "annot_score", "function"][:len(u.columns)]
u = u.dropna(subset=["gene"]).drop_duplicates("gene")
u["annot_score"] = pd.to_numeric(u.annot_score, errors="coerce")
u["no_func"] = u["function"].isna().astype(int)
sel = u[u.gene.isin(d.loc[prom_mask, "gene"])]
ctl = u[u.gene.isin(bg.gene)]
if len(sel) >= 5 and len(ctl) >= 20:
    _, p_score = mannwhitneyu(sel.annot_score.dropna(), ctl.annot_score.dropna())
    or_f, p_func = fisher_exact(
        [[int(sel.no_func.sum()), len(sel) - int(sel.no_func.sum())],
         [int(ctl.no_func.sum()), len(ctl) - int(ctl.no_func.sum())]],
        alternative="greater")
    axes_rows += [
        ("UniProt annotation score (median)",
         "%.1f vs %.1f" % (sel.annot_score.median(), ctl.annot_score.median()),
         None, round(p_score, 4)),
        ("UniProt empty function field",
         "%.1f%% vs %.1f%%" % (100 * sel.no_func.mean(), 100 * ctl.no_func.mean()),
         round(or_f, 3), round(p_func, 4)),
    ]
else:
    print("   WARNING: too few UniProt matches -- the two non-GO dark-proteome "
          "axes are omitted")

pd.DataFrame(
    [(a, c, o, p, "significant" if p is not None and p < 0.05
      else "not significant") for a, c, o, p in axes_rows],
    columns=["axis", "comparison", "odds_ratio", "p", "verdict"]).assign(
    background="coverage-matched non-discordant (n=%d)" % len(bg)).to_csv(
    os.path.join(FIGDATA, "discordance_dark_proteome_axes.csv"), index=False)

print("dark proteome: %d/%d promoted (%.1f%%) carry no BP term"
      % (k_d, n_d, 100 * k_d / n_d))
print("wrote 8 cache tables into figure_data/")
print("   now re-run `make tables` to fold them into discordance_stats.json")
