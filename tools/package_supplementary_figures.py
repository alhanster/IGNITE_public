#!/usr/bin/env python
"""Packages the three supplementary display figures from final_plots/.

Invoked by `make final-outputs`, writing to final_outputs/supplementary_figures/ under the names
it publishes, and by `make supplementary`, writing to supplementary/figures/ under
"Supp Figure <N> - <title>.png". tools/supplementary_figures.tsv is the single source for
the S-numbers.

Usage: .venv/bin/python tools/package_supplementary_figures.py [--out-dir DIR] [--naming STYLE]
"""
import argparse
import os
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
assert os.path.exists(os.path.join(ROOT, "Makefile")), f"not a repo root: {ROOT}"

MANIFEST = os.path.join(ROOT, "tools", "supplementary_figures.tsv")
PLOTS = os.path.join(ROOT, "final_plots")

ILLEGAL = set('/\\:*?"<>|')


def fail(msg):
    raise SystemExit(f"package_supplementary_figures: {msg}")


def read_manifest():
    """Parses the manifest; shares its shape and failure modes with the tables packager's parser."""
    rows = []
    with open(MANIFEST, newline="") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 4:
                fail(f"{MANIFEST}:{lineno}: expected 4 tab-separated fields, got {len(parts)}. "
                     f"Are the separators tabs and not spaces?")
            number, source, title, submission_name = parts[:4]
            if not number.isdigit():
                fail(f"{MANIFEST}:{lineno}: number must be a positive integer, got {number!r}")
            rows.append({"number": int(number), "source": source, "title": title,
                         "submission_name": submission_name, "lineno": lineno})
    return rows


def validate(rows):
    if not rows:
        fail(f"{MANIFEST} lists no figures")

    seen = {}
    for r in rows:
        if r["number"] in seen:
            fail(f"{MANIFEST}:{r['lineno']}: figure number {r['number']} is already used on "
                 f"line {seen[r['number']]}")
        seen[r["number"]] = r["lineno"]

    # A numbering gap means a figure was dropped without renumbering, mismatching the manuscript's cross-references.
    if sorted(seen) != list(range(1, len(rows) + 1)):
        fail(f"figure numbers must run 1..{len(rows)} with no gaps; got {sorted(seen)}")

    for r in rows:
        for field in ("title", "submission_name"):
            if not r[field]:
                fail(f"{MANIFEST}:{r['lineno']}: {field} is empty")
        bad = ILLEGAL & set(r["title"])
        if bad:
            fail(f"{MANIFEST}:{r['lineno']}: title contains {sorted(bad)}, illegal in a filename")
        src = os.path.join(PLOTS, r["source"])
        if not os.path.isfile(src):
            fail(f"{MANIFEST}:{r['lineno']}: no such file: final_plots/{r['source']}"
                 f"\n  Run `make figures` (or `make supplementary`) first -- these are rendered.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join("final_outputs", "supplementary_figures"),
                    help="destination, relative to the repo root")
    ap.add_argument("--naming", choices=["submission", "numbered"], default="submission",
                    help="'submission' keeps the published S<N>_ names; 'numbered' writes "
                         "'Supp Figure <N> - <title>.png'")
    args = ap.parse_args()

    rows = read_manifest()
    validate(rows)

    outd = os.path.join(ROOT, args.out_dir)
    # Output directory is rebuilt from scratch each run so stale filenames from a prior run do not linger.
    if os.path.isdir(outd):
        shutil.rmtree(outd)
    os.makedirs(outd)

    names = {}
    for r in sorted(rows, key=lambda r: r["number"]):
        if args.naming == "submission":
            name = r["submission_name"]
        else:
            ext = os.path.splitext(r["source"])[1]
            name = f"Supp Figure {r['number']} - {r['title']}{ext}"
        names[r["number"]] = name
        shutil.copyfile(os.path.join(PLOTS, r["source"]), os.path.join(outd, name))
        print(f"  Supp Figure {r['number']}: {r['source']} -> {name}")

    print(f"wrote {len(rows)} supplementary figures to {args.out_dir}/")
    copy_pdfs(rows, names, outd, args.out_dir)


def pdf_source(row):
    """PDF path derives from the PNG path (same stem, final_plots/pdf/) to avoid a second source of truth."""
    stem = os.path.splitext(os.path.basename(row["source"]))[0]
    return os.path.join(PLOTS, "pdf", f"{stem}.pdf")


def copy_pdfs(rows, names, outd, out_dir_label):
    """Packages the vector twins into <out-dir>/pdf/, named to match their PNGs.

    See REPRODUCIBILITY.md, PDF outputs and the supplementary target.
    """
    present = [r for r in rows if os.path.isfile(pdf_source(r))]
    if not present:
        print(f"  no PDFs in final_plots/pdf/ -- skipping {out_dir_label}/pdf/."
              f"\n  Expected off macOS: save_figure() writes the PNG, then stops rather than"
              f"\n  downgrading the PDF device. See src/figures/vector_output.R.")
        return
    if len(present) != len(rows):
        missing = [os.path.relpath(pdf_source(r), ROOT) for r in rows if r not in present]
        fail("final_plots/pdf/ holds some but not all supplementary PDFs, so a render is "
             "incomplete rather than absent:\n  missing: " + ", ".join(missing) +
             "\n  Re-run `make figures` (or `make supplementary`) before packaging.")

    pdfd = os.path.join(outd, "pdf")
    os.makedirs(pdfd)
    for r in sorted(rows, key=lambda r: r["number"]):
        name = f"{os.path.splitext(names[r['number']])[0]}.pdf"
        shutil.copyfile(pdf_source(r), os.path.join(pdfd, name))
        print(f"  Supp Figure {r['number']}: pdf/{os.path.basename(pdf_source(r))} -> pdf/{name}")

    print(f"wrote {len(rows)} supplementary PDFs to {out_dir_label}/pdf/")


if __name__ == "__main__":
    main()
