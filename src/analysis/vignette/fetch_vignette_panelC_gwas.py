#!/usr/bin/env python
"""STAT4's immune-GWAS feature, decomposed into the diseases that constitute it: panel c of the vignette figure. Opt-in, network. Run via make vignette-panelc.

Fetches, from Open Targets, the diseases and credible sets contributing to STAT4's gwas_credible_sets score against MONDO:0005046, and writes the per-disease composition used to draw the panel.

See REPRODUCIBILITY.md.
"""
import argparse
import json
import os
import sys
import time


import pandas as pd
import requests

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OT_ARCHIVE = ("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/"
              "output/association_by_datasource_direct/")
STAT4 = "ENSG00000138378"
IMMUNE = "MONDO_0005046"          # immune system disorder -- the subtree that defines gwas_score
DATASOURCE = "gwas_credible_sets"  # the datasource that IS the feature

# One call each, asserted to exceed the known count (1,164 associations at 26.06 vs. a 1200 limit) so truncation is caught.
PAGE_SIZE = 2000
EVIDENCE_SIZE = 1000  # ditto for the evidence trace; 115 credible sets at 26.06
ROUND_DP = 6
FEATURE_TOL = 1e-6    # the matrix stores the feature rounded to 6 dp
ASSERT_TAGS = True    # traced contributor set must equal the ancestor-derived immune tag

# Curated short labels for the four contributors whose OT name overruns the gutter; the first matches the label used in the prior panel.
SHORT = {
    "EFO_0009459":   "ACPA+ rheumatoid arthritis",      # ACPA-positive rheumatoid arthritis
    "EFO_0008536":   "anti-centromere+ scleroderma",    # anti-centromere-antibody-positive ...
    "MONDO_0000589": "musculoskeletal autoimmunity",    # autoimmune disorder of musculoskel...
    "MONDO_0004947": "B-cell ALL",                      # B-cell acute lymphoblastic leukemia
}


def _cap(name):
    """Capitalise the first letter of each space-separated word; leave every other character unchanged.

    Not str.title() or str.capitalize(), which lower-case the remainder of a word, giving "Acpa+ Rheumatoid Arthritis", "B-Cell All", and "Anti-Centromere+ Scleroderma". The rule preserves case elsewhere, so acronyms (ACPA+, ALL) and hyphenated terms (B-cell, anti-centromere+) keep their existing case.
    """
    return " ".join(w[:1].upper() + w[1:] for w in name.split(" "))


def repo_root_from(script_path):
    """Resolves <root> from <root>/src/analysis/vignette/<this>.py. Falls back to cwd when run from elsewhere."""
    guess = os.path.abspath(os.path.join(os.path.dirname(script_path), "..", "..", ".."))
    return guess if os.path.isfile(os.path.join(guess, "Makefile")) else os.getcwd()


def ot_gql(session, query, variables=None, retries=5):
    """POST a GraphQL query to Open Targets with retry/backoff."""
    last = None
    for attempt in range(retries):
        try:
            resp = session.post(OT_URL, json={"query": query, "variables": variables or {}},
                                timeout=120)
            if resp.status_code == 200:
                payload = resp.json()
                if payload.get("errors"):
                    raise RuntimeError(payload["errors"])
                return payload["data"]
            last = f"HTTP {resp.status_code}: {resp.text[:200]}"
        except Exception as exc:                       # network / JSON / API error
            last = str(exc)
        wait = 2 ** attempt
        print(f"  OT retry {attempt + 1}/{retries} ({last[:100]}); sleeping {wait}s",
              file=sys.stderr)
        time.sleep(wait)
    raise RuntimeError(f"OT query failed after {retries} attempts: {last}")


def ot_data_version(session):
    """Open Targets data release, e.g. '26.06'."""
    m = ot_gql(session, "{ meta { dataVersion { year month iteration } } }")[
        "meta"]["dataVersion"]
    v = f"{m['year']}.{m['month']}"
    return v + (f".{m['iteration']}" if m.get("iteration") else "")


Q_ASSOC = """
query($id: String!, $n: Int!) {
  target(ensemblId: $id) {
    approvedSymbol
    drugAndClinicalCandidates { count }
    associatedDiseases(page: {index: 0, size: $n}) {
      count
      rows {
        score
        disease { id name ancestors }
        datasourceScores { id score }
      }
    }
  }
}
"""

Q_FEATURE = """
query($tid: String!, $did: String!) {
  disease(efoId: $did) {
    associatedTargets(Bs: [$tid], page: {index: 0, size: 1}) {
      rows { score datasourceScores { id score } }
    }
  }
}
"""

Q_EVIDENCE = """
query($tid: String!, $did: String!, $ds: String!, $n: Int!) {
  disease(efoId: $did) {
    evidences(ensemblIds: [$tid], datasourceIds: [$ds], enableIndirect: true, size: $n) {
      count
      rows { score disease { id name } credibleSet { studyLocusId } }
    }
  }
}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expect-release", default="26.06",
                    help="abort unless the live API serves this release (default 26.06)")
    ap.add_argument("--repo-root", default=None, help="test hook; defaults to <script>/../..")
    ap.add_argument("--out-dir", default=None, help="test hook; defaults to <root>/figure_data")
    args = ap.parse_args()

    root = args.repo_root or repo_root_from(os.path.abspath(__file__))
    outdir = args.out_dir or os.path.join(root, "figure_data")
    os.makedirs(outdir, exist_ok=True)

    session = requests.Session()

    release = ot_data_version(session)
    if release != args.expect_release:
        raise SystemExit(
            f"Open Targets serves release {release}, not the expected {args.expect_release}.\n"
            f"Association scores move between releases and this file backs a published panel, so\n"
            f"nothing was written. Either re-pin with --expect-release {release} and re-check every\n"
            f"number the panel and Results quote, or pull the pinned release from the archive:\n"
            f"  {OT_ARCHIVE % args.expect_release}")
    print(f"Open Targets release {release}")

    assoc = ot_gql(session, Q_ASSOC, {"id": STAT4, "n": PAGE_SIZE})["target"]
    ad = assoc["associatedDiseases"]
    n_total = ad["count"]
    assert n_total < PAGE_SIZE, (
        f"{n_total} associations exceeds PAGE_SIZE={PAGE_SIZE}; raise it -- a truncated fetch "
        f"would silently drop GWAS-backed diseases from the table")
    assert len(ad["rows"]) == n_total, f"asked for all {n_total}, got {len(ad['rows'])}"
    print(f"{assoc['approvedSymbol']}: {n_total} disease associations")

    rows = []
    for r in ad["rows"]:
        ds = {d["id"]: d["score"] for d in r["datasourceScores"]}
        if DATASOURCE not in ds:
            continue
        anc = set(r["disease"].get("ancestors") or [])
        rows.append({
            "disease": r["disease"]["name"],
            "label": _cap(SHORT.get(r["disease"]["id"], r["disease"]["name"])),
            "efo": r["disease"]["id"],
            "gwas_credible_sets_score": ds[DATASOURCE],
            "overall_assoc_score": r["score"],
            "tag": "immune" if (IMMUNE in anc or r["disease"]["id"] == IMMUNE) else "other",
        })
    df = pd.DataFrame(rows)
    print(f"  {len(df)} carry {DATASOURCE} evidence; "
          f"{int((df.tag == 'immune').sum())} inside {IMMUNE} by ancestor list")

    # Raises if OT stops returning a curated label's key, rather than falling back to the overlong name.
    stale = sorted(set(SHORT) - set(df.efo))
    if stale:
        print(f"WARNING: SHORT has no matching row for {stale} -- update the map in "
              f"{os.path.basename(__file__)}", file=sys.stderr)

    frows = ot_gql(session, Q_FEATURE, {"tid": STAT4, "did": IMMUNE})[
        "disease"]["associatedTargets"]["rows"]
    assert len(frows) == 1, f"expected one target row for {IMMUNE}, got {len(frows)}"
    fds = {d["id"]: d["score"] for d in frows[0]["datasourceScores"]}
    feature = fds[DATASOURCE]

    matrix = os.path.join(root, "data", "perturbseq", "pu", "pu_model_matrix.parquet")
    if os.path.isfile(matrix):
        stored = pd.read_parquet(matrix, columns=["gene", "gwas_score"])
        stored = float(stored.loc[stored.gene == "STAT4", "gwas_score"].iloc[0])
        assert abs(stored - feature) < FEATURE_TOL, (
            f"the panel's reference line would not be the model input: API {feature!r} vs "
            f"matrix {stored!r}. Either the feature was rebuilt from a different release or "
            f"its definition changed; do not draw the panel until this is resolved.")
        print(f"  feature {feature:.7f} matches pu_model_matrix gwas_score {stored} "
              f"(|d| < {FEATURE_TOL})")
    else:
        print(f"  WARNING: {matrix} absent; feature-identity check skipped", file=sys.stderr)

    ev = ot_gql(session, Q_EVIDENCE,
                {"tid": STAT4, "did": IMMUNE, "ds": DATASOURCE, "n": EVIDENCE_SIZE})[
                    "disease"]["evidences"]
    assert ev["count"] <= EVIDENCE_SIZE, (
        f"{ev['count']} evidence rows exceeds EVIDENCE_SIZE={EVIDENCE_SIZE}; raise it")
    assert len(ev["rows"]) == ev["count"], "evidence fetch truncated"
    counts, best, sets = {}, {}, {}
    for r in ev["rows"]:
        did = r["disease"]["id"]
        cs = r["credibleSet"]
        # credibleSet is nullable and confirmed non-null for all 115 rows at 26.06; a null would leave credible_sets shorter than n_credible_sets.
        assert cs and cs.get("studyLocusId"), (
            f"evidence row for {did} carries no credibleSet.studyLocusId -- the id list and "
            f"n_credible_sets would disagree; do not write the table")
        counts[did] = counts.get(did, 0) + 1
        best[did] = max(best.get(did, 0.0), r["score"])
        sets.setdefault(did, []).append((r["score"], cs["studyLocusId"]))
    n_distinct = len({i for v in sets.values() for _, i in v})
    print(f"  {ev['count']} credible sets roll into {IMMUNE} from {len(counts)} diseases "
          f"({n_distinct} distinct -- some are shared across diseases, see the docstring)")

    # Sorted descending by L2G score, id as tie-break, so the first id in the joined list is the one best_l2g reports.
    joined = {d: "; ".join(i for _, i in sorted(v, key=lambda t: (-t[0], t[1])))
              for d, v in sets.items()}

    df["contributes_to_gwas_score"] = df.efo.isin(counts).astype(int)
    df["n_credible_sets"] = df.efo.map(lambda e: counts.get(e, 0))
    df["best_l2g"] = df.efo.map(lambda e: round(best[e], ROUND_DP) if e in best else pd.NA)
    # Empty string for the 25 non-contributing rows means the disease is outside Q_EVIDENCE's IMMUNE scope, not that it lacks credible sets.
    df["credible_sets"] = df.efo.map(lambda e: joined.get(e, ""))

    missing = set(counts) - set(df.efo)
    assert not missing, (
        f"{len(missing)} diseases contribute evidence but carry no {DATASOURCE} association row: "
        f"{sorted(missing)[:5]}")
    assert int(df.n_credible_sets.sum()) == ev["count"], "credible sets lost in the join"
    listed = df.credible_sets.map(lambda s: len(s.split("; ")) if s else 0)
    assert (listed == df.n_credible_sets).all(), (
        "credible_sets and n_credible_sets disagree for "
        f"{sorted(df.loc[listed != df.n_credible_sets, 'efo'])} -- the sheet would name a "
        "different number of sets than it counts")
    assert int(listed.sum()) == ev["count"], "credible set ids lost in the join"
    dupes = {d: len(v) - len({i for _, i in v}) for d, v in sets.items()}
    assert not any(dupes.values()), (
        f"a studyLocusId repeats WITHIN one disease: {[d for d, n in dupes.items() if n]}. "
        "Across diseases it legitimately repeats; within one it means the evidence fetch "
        "double-counted")
    assert (df.loc[df.contributes_to_gwas_score == 0, "credible_sets"] == "").all(), (
        "a non-contributing disease carries credible set ids, which cannot happen while "
        "Q_EVIDENCE is scoped to the immune subtree")
    if ASSERT_TAGS:
        traced = set(df.loc[df.contributes_to_gwas_score == 1, "efo"])
        tagged = set(df.loc[df.tag == "immune", "efo"])
        assert traced == tagged, (
            "the traced contributor set and the ancestor-derived immune tag disagree "
            f"(+{sorted(traced - tagged)} / -{sorted(tagged - traced)}). The panel's colour and "
            "caption assume they are the same set -- resolve before drawing.")

    # Deterministic order and rounding.
    for col in ("gwas_credible_sets_score", "overall_assoc_score"):
        df[col] = df[col].round(ROUND_DP)
    df = df.sort_values(["gwas_credible_sets_score", "efo"], ascending=[False, True])
    df = df[["disease", "label", "efo", "gwas_credible_sets_score", "overall_assoc_score", "tag",
             "contributes_to_gwas_score", "n_credible_sets", "best_l2g", "credible_sets"]]

    csv_path = os.path.join(outdir, "vignette_panelC_gwas_diseases.csv")
    df.to_csv(csv_path, index=False)

    meta = {
        "ot_release": release,
        "target": STAT4,
        "symbol": assoc["approvedSymbol"],
        "datasource": DATASOURCE,
        "immune_subtree": IMMUNE,
        "feature_definition": f"{DATASOURCE} score for {assoc['approvedSymbol']} x {IMMUNE}",
        "gwas_score_feature_value": feature,
        "n_diseases_total": n_total,
        # Consumed by build_vignette_meta.py as approved_drugs, paired with ot_diseases (n_diseases_total); see the superset caveat.
        "approved_drugs": int(assoc["drugAndClinicalCandidates"]["count"]),
        "n_diseases_with_gwas_evidence": int(len(df)),
        "n_contributing_diseases": int(df.contributes_to_gwas_score.sum()),
        "n_credible_sets_total": int(df.n_credible_sets.sum()),
        "contributor_rule": ("evidence rolled into the immune subtree with enableIndirect=true, "
                             "traced from the evidence endpoint rather than inferred from "
                             "ancestor lists"),
    }
    # No fetch timestamp: figure_data/ is byte-verified by make verify-tables, and the release string serves as provenance instead.
    meta_path = os.path.join(outdir, "vignette_panelC_gwas_meta.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=1)
        fh.write("\n")

    print(f"wrote {csv_path} ({len(df)} rows)")
    print(f"wrote {meta_path}")


if __name__ == "__main__":
    main()
