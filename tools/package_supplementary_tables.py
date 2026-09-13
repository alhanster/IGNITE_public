#!/usr/bin/env python
"""Package the numbered supplementary tables from figure_data/.

Called by make final-outputs (writes to final_outputs/supplementary_tables/, the default) and by make supplementary (writes to supplementary/tables/, via --out-dir). Reads tools/supplementary_tables.tsv, the manifest that sets each table's number and display name, and copies each source file under its manuscript name. The folder README only points to Supplementary Information.pdf, which describes each table.

See REPRODUCIBILITY.md, Numbered supplementary tables, for the design rationale.

Run: PYTHONPATH=src .venv/bin/python tools/package_supplementary_tables.py
"""
import argparse
import os
import shutil
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"

MANIFEST = os.path.join(ROOT, "tools", "supplementary_tables.tsv")

# Also in README.md, REPRODUCIBILITY.md and tools/write_submission_readme.py; keep the four in step.
PLATFORM_STATEMENT = (
    'All results reported in the manuscript and supplementary tables were generated on macOS '
    '(Apple Silicon) using the environment provided in this repository. Minor numerical '
    'differences may occur when the analyses are reproduced on other operating systems or '
    'hardware, because of differences in numerical libraries, floating-point operations and '
    'parallel computation. These differences do not affect the qualitative conclusions or '
    'statistical interpretation of the analyses. In the exploratory discordance analysis '
    '(Supplementary Fig. S2 and Supplementary Table 9), the core gene set and features near '
    'the FDR threshold may differ.')

# Sourced from figure_data/, not final_plots/supplementary/. See REPRODUCIBILITY.md, Numbered supplementary tables, for rationale.
SUPP = os.path.join(ROOT, "figure_data")
# Default directory preserves make final-outputs' prior behavior; make supplementary overrides it via --out-dir.
DEFAULT_OUT = os.path.join("final_outputs", "supplementary_tables")

# A supplementary source with no producer script is never assigned a table number. See PROVENANCE.md, Written by stage 1 (87 files), for what each table is.

# Filenames with characters illegal for submission portals or Windows are flagged as an error, not rewritten automatically.
ILLEGAL = set('/\\:*?"<>|')


def fail(msg):
    raise SystemExit(f"package_supplementary_tables: {msg}")


def read_manifest():
    rows = []
    with open(MANIFEST, newline="") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                fail(f"{MANIFEST}:{lineno}: expected 4 tab-separated fields, got {len(parts)}. "
                     f"Are the separators tabs and not spaces?")
            number, source, title, description = (p.strip() for p in parts[:4])
            if not number.isdigit():
                fail(f"{MANIFEST}:{lineno}: number must be a positive integer, got {number!r}")
            rows.append({"number": int(number), "title": title, "description": description,
                         "lineno": lineno, **parse_source(source, lineno)})
    return rows


def parse_source(source, lineno):
    """One file, or a workbook of sheets combined at packaging time; excluded from figure_data/ and final_plots/ because zip timestamps make xlsx non-byte-reproducible."""
    if "=" not in source:
        return {"source": source, "sheets": None}
    sheets = {}
    for part in source.split(";"):
        if "=" not in part:
            fail(f"{MANIFEST}:{lineno}: sheet spec {part!r} is not 'Sheet name=file.csv'")
        name, fname = (x.strip() for x in part.split("=", 1))
        if not name or not fname:
            fail(f"{MANIFEST}:{lineno}: sheet spec {part!r} has an empty name or filename")
        # Excel's own format limits on size and name, checked here so a failure clearly names Excel as the cause.
        if len(name) > 31:
            fail(f"{MANIFEST}:{lineno}: sheet name {name!r} exceeds Excel's 31-character limit")
        if set(name) & set('[]:*?/\\'):
            fail(f"{MANIFEST}:{lineno}: sheet name {name!r} contains a character Excel forbids")
        if name in sheets:
            fail(f"{MANIFEST}:{lineno}: duplicate sheet name {name!r}")
        sheets[name] = fname
    return {"source": None, "sheets": sheets}


def validate(rows):
    if not rows:
        fail(f"{MANIFEST} lists no tables")

    seen = {}
    for r in rows:
        if r["number"] in seen:
            fail(f"{MANIFEST}:{r['lineno']}: table number {r['number']} is already used on "
                 f"line {seen[r['number']]}")
        seen[r["number"]] = r["lineno"]

    # A numbering gap means a table was dropped without renumbering, which would mismatch the manuscript's cross-references; raised as an error, not passed through silently.
    expected = list(range(1, len(rows) + 1))
    if sorted(seen) != expected:
        fail(f"table numbers must run 1..{len(rows)} with no gaps; got {sorted(seen)}")

    for r in rows:
        if not r["title"]:
            fail(f"{MANIFEST}:{r['lineno']}: title is empty")
        bad = ILLEGAL & set(r["title"])
        if bad:
            fail(f"{MANIFEST}:{r['lineno']}: title contains {sorted(bad)}, illegal in a filename")
        if not r["description"]:
            fail(f"{MANIFEST}:{r['lineno']}: description is empty -- every numbered table needs one")

        for fname in sources_of(r):
            src = os.path.join(SUPP, fname)
            if not os.path.isfile(src):
                fail(f"{MANIFEST}:{r['lineno']}: no such file: figure_data/{fname}"
                     f"\n  Run `make tables` first -- these are stage-1 outputs.")


def sources_of(row):
    """Every file a manifest row depends on, for both single-file and multi-sheet rows."""
    return [row["source"]] if row["sheets"] is None else list(row["sheets"].values())


def write_workbook(sheets, out_path):
    """Combines per-sheet CSVs into one xlsx, sheet order set by the manifest; the only place an xlsx is produced, written outside every hashed tree."""
    import pandas as pd  # local: only the multi-sheet path needs it

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for name, path in sheets.items():
            # dtype=str preserves the CSV's formatting, so blank CI bounds on enrichment reference rows stay blank and no numeric value is re-rendered.
            pd.read_csv(path, dtype=str, keep_default_na=False).to_excel(
                writer, sheet_name=name, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=DEFAULT_OUT,
                    help="destination, relative to the repo root")
    args = ap.parse_args()
    OUTD = os.path.join(ROOT, args.out_dir)

    rows = read_manifest()
    validate(rows)

    # Output directory is deleted and rebuilt each run so renumbering cannot leave stale filenames behind, as make final-outputs does for final_outputs/.
    if os.path.isdir(OUTD):
        shutil.rmtree(OUTD)
    os.makedirs(OUTD)

    for r in sorted(rows, key=lambda r: r["number"]):
        multi = r["sheets"] is not None
        ext = ".xlsx" if multi else os.path.splitext(r["source"])[1]
        name = f"Supp Table {r['number']} - {r['title']}{ext}"
        dest = os.path.join(OUTD, name)

        if multi:
            write_workbook({s: os.path.join(SUPP, f) for s, f in r["sheets"].items()}, dest)
        else:
            shutil.copyfile(os.path.join(SUPP, r["source"]), dest)

        print(f"  Supp Table {r['number']}: -> {name}")

    with open(os.path.join(OUTD, "README.md"), "w") as fh:
        fh.write("# Supplementary tables\n\n"
                 "See Supplementary Information.pdf for a description of each table in this folder.\n")

    print(f"wrote {len(rows)} numbered tables to {args.out_dir}/")


if __name__ == "__main__":
    main()
