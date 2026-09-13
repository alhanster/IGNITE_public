#!/usr/bin/env python3
"""Shared gene-name utilities: formatting checks and HGNC harmonization."""

from __future__ import annotations

import re
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/tsv/"
            "tsv/hgnc_complete_set.txt")
HGNC_CACHE = Path(tempfile.gettempdir()) / "hgnc_complete_set.txt"

# The pinned reference the build harmonizes against. Named hgnc_symbol_subset.txt here (a
# four-column projection) and hgnc_complete_set.txt in the upstream feature-build repo.
_SOURCES = Path(__file__).resolve().parents[3] / "data" / "perturbseq" / "sources"
PINNED_HGNC = next((p for p in (_SOURCES / "hgnc_symbol_subset.txt",
                                _SOURCES / "hgnc_complete_set.txt") if p.exists()), None)

# Body may be lowercase for orf names (C1orf112) and hyphenated readthrough/antisense symbols (ANKRD13C-DT).
FORMAT_RE = re.compile(r"^[A-Z0-9][A-Za-z0-9._-]*$")


# ---------------------------------------------------------------------------
def load_hgnc(path: str | Path | None = None, refresh: bool = False,
              allow_download: bool = False) -> dict:
    """Load the pinned HGNC reference and derive the rename maps from it.

    Returns {"approved": set, "prev2app": {prev_symbol: approved}, "alias": set,
    "n_records": int, "release": str}, built from records with status Approved.

    A missing local reference raises rather than falling back to the network. HGNC renames
    symbols between releases, so a silent fetch harmonizes to a different symbol space than
    the rest of the committed data and re-keys every downstream join on `gene`. Pass
    allow_download=True only to deliberately refresh against the current release.
    """
    if path is not None:
        src = Path(path)
        if not src.exists():
            raise FileNotFoundError(f"HGNC file not found: {src}")
    elif not allow_download:
        raise FileNotFoundError(
            "load_hgnc() needs an explicit path to the pinned HGNC reference "
            "(data/perturbseq/sources/hgnc_symbol_subset.txt). Refusing to fetch the current "
            "release: HGNC renames symbols between releases, so downloading here would "
            "harmonize to a different symbol space than the committed data and silently "
            "re-key every gene join. Pass allow_download=True to override deliberately.")
    else:
        src = HGNC_CACHE
        if refresh or not src.exists():
            print(f"[hgnc] downloading complete set -> {src}")
            urllib.request.urlretrieve(HGNC_URL, src)
        else:
            print(f"[hgnc] using cached set {src}")

    # date_modified is absent from the committed column subset and present in the upstream
    # complete set; it stamps the release when available.
    need = ["symbol", "status", "alias_symbol", "prev_symbol"]
    have = pd.read_csv(src, sep="\t", nrows=0).columns
    h = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False,
                    usecols=need + [c for c in ["date_modified"] if c in have])
    h = h[h["status"] == "Approved"]

    dates = h["date_modified"][h["date_modified"] != ""] if "date_modified" in h else None
    release = f"latest date_modified {dates.max()}" if dates is not None and len(dates) \
        else "date_modified not carried by this file"
    print(f"[hgnc] {src}: {len(h)} approved records, {release}")

    approved = set(h["symbol"])
    prev2app: dict[str, str] = {}
    alias: set[str] = set()
    for sym, prev, ali in zip(h["symbol"], h["prev_symbol"], h["alias_symbol"]):
        for p in prev.split("|"):
            if p:
                prev2app.setdefault(p, sym)
        for a in ali.split("|"):
            if a:
                alias.add(a)
    return {"approved": approved, "prev2app": prev2app, "alias": alias,
            "n_records": len(h), "release": release}


# ---------------------------------------------------------------------------
def check_format(genes, verbose: bool = True) -> dict:
    """Gene-name formatting report, no network or HGNC lookup involved."""
    g = pd.Series(list(genes), dtype="object").astype(str)
    nonstd = sorted(g[~g.str.match(FORMAT_RE)].unique())
    report = {
        "n": len(g),
        "unique": int(g.nunique()),
        "whitespace": int((g != g.str.strip()).sum()),
        "space": int(g.str.contains(" ").sum()),
        "empty_or_na": int(((g.str.strip() == "") | (g == "NA")).sum()),
        "duplicates": int(g.duplicated().sum()),
        "nonstandard": nonstd,
        "special_chars": sorted({c for tok in g for c in tok if not c.isalnum()}),
    }
    if verbose:
        print("--- gene-name formatting check ---")
        print(f"genes: {report['n']} | unique: {report['unique']}")
        print(f"leading/trailing whitespace: {report['whitespace']}")
        print(f"embedded space: {report['space']}")
        print(f"empty or 'NA': {report['empty_or_na']}")
        print(f"duplicate symbols: {report['duplicates']}")
        print(f"non-standard-shape symbols: {len(report['nonstandard'])}"
              + (f" -> {report['nonstandard'][:15]}" if report["nonstandard"] else ""))
        print(f"special characters present: {report['special_chars']}")
    return report


# ---------------------------------------------------------------------------
def classify_hgnc(genes, hgnc: dict, verbose: bool = True) -> pd.Series:
    """Classify each gene as approved, previous_symbol, alias_symbol, or not_found."""
    approved, prev2app, alias = hgnc["approved"], hgnc["prev2app"], hgnc["alias"]

    def cl(x: str) -> str:
        if x in approved:
            return "approved"
        if x in prev2app:
            return "previous_symbol"
        if x in alias:
            return "alias_symbol"
        return "not_found"

    g = pd.Series(list(genes), dtype="object").astype(str)
    st = g.map(cl)
    if verbose:
        print("--- HGNC status ---")
        vc = st.value_counts()
        for k in ["approved", "previous_symbol", "alias_symbol", "not_found"]:
            if k in vc:
                print(f"{k}: {int(vc[k])}")
        print(f"approved %: {round(100 * (st == 'approved').mean(), 1)}")
    return st


# ---------------------------------------------------------------------------
def harmonize(df: pd.DataFrame, hgnc: dict, gene_col: str = "gene",
              verbose: bool = True):
    """Renames previous symbols to approved ones, dropping colliding rows; info has renamed count, dropped list."""
    prev2app, approved = hgnc["prev2app"], hgnc["approved"]
    present = set(df[gene_col])

    def target(x: str) -> str:
        return prev2app[x] if (x not in approved and x in prev2app) else x

    tgt = df[gene_col].map(target)
    is_prev = tgt != df[gene_col]
    collide = is_prev & tgt.isin(present)          # target already exists → drop
    dropped = sorted(df.loc[collide, gene_col].unique())

    out = df.loc[~collide].copy()
    out[gene_col] = out[gene_col].map(target)
    n_renamed = int((is_prev & ~collide).sum())

    # Guard: two previous symbols mapping to the same new target.
    if out[gene_col].duplicated().any():
        dups = sorted(out.loc[out[gene_col].duplicated(keep=False), gene_col].unique())
        raise ValueError(f"harmonize produced duplicate gene(s): {dups}")

    info = {"renamed": n_renamed, "dropped": dropped}
    if verbose:
        print("--- HGNC harmonization ---")
        print(f"symbols renamed: {n_renamed} | deprecated dropped: {len(dropped)}"
              + (f" ({', '.join(dropped)})" if dropped else ""))
    return out, info


# ---------------------------------------------------------------------------
def _main(argv: list[str]) -> None:
    if not argv:
        sys.exit("usage: python src/analysis/common/_gene_names.py <csv> [gene_col] [--download]")
    download = "--download" in argv          # opt in to the current release, off by default
    argv = [a for a in argv if a != "--download"]
    path = argv[0]
    gene_col = argv[1] if len(argv) > 1 else "gene"
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    if gene_col not in df.columns:
        sys.exit(f"column '{gene_col}' not in {path}; columns: {list(df.columns)}")
    genes = df[gene_col]
    check_format(genes)
    print()
    classify_hgnc(genes, load_hgnc(None if download else PINNED_HGNC,
                                   allow_download=download))


if __name__ == "__main__":
    _main(sys.argv[1:])
