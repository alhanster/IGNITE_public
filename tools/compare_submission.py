#!/usr/bin/env python
"""Compares two assembled final_outputs/ trees against a reference package, ignoring fields
that cannot match. See REPRODUCIBILITY.md.

python3 tools/compare_submission.py <reference-package-dir> [final_outputs]
"""
import hashlib, os, re, sys, zipfile

XLSX_IGNORED_PARTS = {"docProps/core.xml"}
PDF_VOLATILE = re.compile(rb"/(CreationDate|ModDate)\s*\(D:[^)]*\)")
PDF_ID = re.compile(rb"/ID\s*\[\s*<[0-9a-fA-F]*>\s*<[0-9a-fA-F]*>\s*\]")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def compare_xlsx(pa, pb):
    """Hash every part of both workbooks; the timestamp part is the only allowed difference."""
    with zipfile.ZipFile(pa) as za, zipfile.ZipFile(pb) as zb:
        na, nb = set(za.namelist()), set(zb.namelist())
        if na != nb:
            return False, f"part list differs: {sorted(na ^ nb)}"
        bad = [p for p in sorted(na)
               if p not in XLSX_IGNORED_PARTS and sha(za.read(p)) != sha(zb.read(p))]
        if bad:
            return False, f"{len(bad)} part(s) differ: {bad}"
        return True, f"{len(na) - len(XLSX_IGNORED_PARTS & na)} parts identical, timestamp only"


def compare_pdf(pa, pb):
    """Blank the date fields and the /ID derived from them, then compare what is left."""
    def norm(p):
        b = open(p, "rb").read()
        return PDF_ID.sub(b"/ID[]", PDF_VOLATILE.sub(b"/DATE()", b))
    a, b = norm(pa), norm(pb)
    if a == b:
        return True, "identical, creation stamp only"
    n = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
    return False, f"{n} byte(s) differ beyond the creation stamp"


# Cell classes exclude \n; otherwise size/hash matches can span rows and blank a real edit as matching.
MANIFEST_ROW = re.compile(
    rb"^(\|\s*`[^`]+\.(?:xlsx|pdf)`\s*\|)[^|\n]*\|[^|\n]*\|", re.MULTILINE)


def compare_manifest(pa, pb):
    """final_outputs/README.md: normalise the rows whose hashes track a creation stamp."""
    def norm(p):
        return MANIFEST_ROW.sub(rb"\1 SIZE | HASH |", open(p, "rb").read())
    a, b = norm(pa), norm(pb)
    if a == b:
        return True, "identical, xlsx/pdf manifest rows only"
    return False, "prose or a png/csv row differs"


def compare_exact(pa, pb):
    a, b = open(pa, "rb").read(), open(pb, "rb").read()
    return (a == b), ("identical (byte-for-byte)" if a == b else
                      f"DIFFERS ({len(a)} vs {len(b)} bytes)")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ref = sys.argv[1]
    cur = sys.argv[2] if len(sys.argv) > 2 else "final_outputs"
    for d in (ref, cur):
        if not os.path.isdir(d):
            sys.exit(f"not a directory: {d}\nRun `make final-outputs` first.")

    def tree(root):
        out = set()
        for dirpath, _, names in os.walk(root):
            for n in names:
                if n == ".DS_Store":
                    continue
                out.add(os.path.relpath(os.path.join(dirpath, n), root))
        return out

    ta, tb = tree(ref), tree(cur)
    ok = True

    for missing, where in ((tb - ta, ref), (ta - tb, cur)):
        for f in sorted(missing):
            print(f"  MISSING from {where}: {f}")
            ok = False

    counts = {}
    for f in sorted(ta & tb):
        ext = os.path.splitext(f)[1].lower()
        pa, pb = os.path.join(ref, f), os.path.join(cur, f)
        if ext == ".xlsx":
            good, note = compare_xlsx(pa, pb)
        elif ext == ".pdf":
            good, note = compare_pdf(pa, pb)
        elif f == "README.md":
            good, note = compare_manifest(pa, pb)
        else:
            good, note = compare_exact(pa, pb)
        counts[ext] = counts.get(ext, [0, 0])
        counts[ext][0 if good else 1] += 1
        if not good:
            ok = False
        print(f"  {'ok  ' if good else 'FAIL'}  {f}  --  {note}")

    print()
    for ext in sorted(counts):
        good, bad = counts[ext]
        print(f"  {ext or '(none)':6s} {good} reproduced" + (f", {bad} DIFFER" if bad else ""))
    print()
    print("SUBMISSION PACKAGE REPRODUCES" if ok else "PACKAGE DIFFERS -- see FAIL lines above")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
