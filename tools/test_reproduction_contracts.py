"""Fast pre-flight checks for the reproduction path.

These do not fit models or render figures. They confirm that the build's contracts hold
before a long run starts: pins parse, the two-stage separation is intact, every file the
submission target copies has a producer, and the manifests are readable.

Run with:  make test
"""
import csv
import hashlib
import os
import re
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src", "analysis", "common"))
sys.path.insert(0, os.path.join(ROOT, "src"))


def _read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return fh.read()


# --- version pins -----------------------------------------------------------------

def test_python_requires_is_parseable():
    import _version_guard as vg
    assert vg._parse_python_requires() is not None, "requirements.txt lost its python_requires directive"


def test_every_declared_pin_is_exact():
    import _version_guard as vg
    pins = vg._parse_exact_pins()
    for pkg in ("numpy", "pandas", "scipy", "scikit-learn", "xgboost", "pyarrow", "matplotlib"):
        assert pkg in pins, f"{pkg} is no longer exactly pinned in requirements.txt"


def test_critical_and_advisory_are_disjoint_and_covered():
    import _version_guard as vg
    assert not set(vg.CRITICAL) & set(vg.ADVISORY)
    pins = set(vg._parse_exact_pins())
    assert pins == set(vg.CRITICAL) | set(vg.ADVISORY), (
        "a pin in requirements.txt is neither enforced nor advisory")


def test_every_xgboost_trainer_is_guarded():
    import _version_guard as vg
    assert vg.audit_coverage() == [], "an XGBClassifier trainer does not call check_pins()"


# --- two-stage contract -----------------------------------------------------------

FIGDIR = os.path.join(ROOT, "src", "figures")


def test_stage2_reads_only_figure_data():
    """Renderers must not reach into data/ or outputs/; stage 2 runs from figure_data/ alone."""
    offenders = []
    for name in sorted(os.listdir(FIGDIR)):
        if not name.endswith(".R"):
            continue
        body = _read(os.path.join("src", "figures", name))
        for bad in ('"data/', "'data/", '"outputs/', "'outputs/"):
            if bad in body:
                offenders.append((name, bad))
    assert offenders == [], f"stage-2 renderer reaches outside figure_data/: {offenders}"


def test_figure_data_tables_referenced_by_renderers_exist():
    missing = []
    for name in sorted(os.listdir(FIGDIR)):
        if not name.endswith(".R"):
            continue
        body = _read(os.path.join("src", "figures", name))
        for ref in set(re.findall(r'["\']([A-Za-z0-9_.\-]+\.(?:csv|json))["\']', body)):
            cand = os.path.join(ROOT, "figure_data", ref)
            if not os.path.exists(cand) and not ref.startswith(("figure_", "Supp")):
                missing.append((name, ref))
    assert missing == [], f"renderer reads a figure_data table that is not committed: {missing}"


# --- submission manifests ---------------------------------------------------------

def _tsv(rel):
    rows = []
    for line in _read(rel).splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.split("\t"))
    return rows


def test_supplementary_manifests_parse():
    """Both manifests must satisfy the contract package_supplementary_*.py enforces:
    at least four tab-separated fields per row, of which the first four are used.
    Trailing empty fields are tolerated, matching the parser's parts[:4] slice."""
    for rel in ("tools/supplementary_figures.tsv", "tools/supplementary_tables.tsv"):
        rows = _tsv(rel)
        assert rows, f"{rel} has no data rows"
        short = [i for i, r in enumerate(rows, 1) if len(r) < 4]
        assert short == [], f"{rel} rows with fewer than 4 fields: {short}"
        blank = [i for i, r in enumerate(rows, 1) if not r[0].strip().isdigit()]
        assert blank == [], f"{rel} rows whose first field is not a number: {blank}"


def test_supplementary_numbers_are_unique_and_contiguous():
    for rel in ("tools/supplementary_figures.tsv", "tools/supplementary_tables.tsv"):
        nums = sorted(int(r[0]) for r in _tsv(rel))
        assert nums == list(range(1, len(nums) + 1)), f"{rel} numbering is not 1..n: {nums}"


# --- build wiring -----------------------------------------------------------------

def test_makefile_declares_the_submission_target():
    mk = _read("Makefile")
    assert re.search(r"^submission:", mk, re.M), "make submission is missing"


def test_every_stage1_script_named_by_the_makefile_exists():
    mk = _read("Makefile")
    named = re.findall(r"src/analysis/((?:[A-Za-z0-9_]+/)?[A-Za-z0-9_]+\.py)", mk)
    assert len(named) >= 30, f"found only {len(named)} script paths; has the path layout changed?"
    missing = [p for p in named if not os.path.exists(os.path.join(ROOT, "src", "analysis", p))]
    assert missing == [], f"Makefile invokes a script that does not exist: {missing}"


def test_every_renderer_named_by_the_makefile_exists():
    mk = _read("Makefile")
    missing = [p for p in re.findall(r"src/figures/([A-Za-z0-9_]+\.R)", mk)
               if not os.path.exists(os.path.join(ROOT, "src", "figures", p))]
    assert missing == [], f"Makefile invokes a renderer that does not exist: {missing}"


def test_no_shell_script_is_syntactically_broken():
    broken = []
    tools = os.path.join(ROOT, "tools")
    for name in sorted(os.listdir(tools)):
        if name.endswith(".sh"):
            r = subprocess.run(["bash", "-n", os.path.join(tools, name)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                broken.append(name)
    assert broken == [], f"shell script fails bash -n: {broken}"


# --- documentation cross-references -----------------------------------------------

# The source tree points at REPRODUCIBILITY.md and figure_data/PROVENANCE.md by document
# name only, never by section, so the documents can be restructured without touching source
# comments. These two tests keep that convention: a pointer that regrows a ", Section" tail
# would couple the two again and drift the moment a heading is renamed.

XREF_DOCS = {"REPRODUCIBILITY.md": "REPRODUCIBILITY.md",
             "PROVENANCE.md": os.path.join("figure_data", "PROVENANCE.md")}

XREF_RE = re.compile(r"(REPRODUCIBILITY|PROVENANCE)\.md,\s*(.{3,160})", re.S)


def _xref_sources():
    out = []
    for sub in ("src", "tools"):
        for dirpath, _, names in os.walk(os.path.join(ROOT, sub)):
            for n in names:
                # This file quotes the pointer pattern in its own comments, so it is not a source.
                if n.endswith((".py", ".R", ".sh", ".md")) and n != os.path.basename(__file__):
                    out.append(os.path.relpath(os.path.join(dirpath, n), ROOT))
    out += ["Makefile", "README.md", ".gitignore", "requirements.txt",
            "src/figures/R-requirements.txt", "REPRODUCIBILITY.md", "data/DATA_AVAILABILITY.md"]
    return [f for f in sorted(set(out)) if os.path.exists(os.path.join(ROOT, f))]


def test_referenced_docs_exist():
    """Every document the tree points at must be present."""
    missing = [rel for rel in sorted(XREF_DOCS.values())
               if not os.path.exists(os.path.join(ROOT, rel))]
    assert missing == [], f"source comments point at absent documents: {missing}"


def test_doc_pointers_do_not_name_sections():
    """A pointer names the document only; a ", Section" tail recouples docs to source."""
    sectioned = []
    for rel in _xref_sources():
        text = _read(rel)
        for m in XREF_RE.finditer(text):
            line = text[:m.start()].count("\n") + 1
            sectioned.append(f"{rel}:{line} -> {m.group(2)[:60]!r}")
    assert sectioned == [], (
        "doc pointer names a section; reference the document only:\n  "
        + "\n  ".join(sectioned))


# --- submission package integrity ---------------------------------------------------

def test_submission_manifest_covers_every_packaged_file():
    """final_outputs/README.md must carry a size and SHA-256 for every file it ships.

    final_outputs/ is gitignored, so a packaged file with no manifest row has no recovery
    path -- the state REPRODUCIBILITY.md's "Files outside the manifest" paragraph says must
    not be reported as protected. Skipped rather than failed on a tree that has not built
    the package, since the rest of this suite runs from a fresh clone.
    """
    sub = os.path.join(ROOT, "final_outputs")
    readme = os.path.join(sub, "README.md")
    if not os.path.isfile(readme):
        pytest.skip("final_outputs/ not built; run `make final-outputs`")

    with open(readme, encoding="utf-8") as fh:
        listed = set(re.findall(r"^\| `([^`]+)` \| [\d,]+ \| `[0-9a-f]{16}` \|$",
                                fh.read(), re.M))
    on_disk = {os.path.relpath(os.path.join(dp, n), sub)
               for dp, _, names in os.walk(sub) for n in names if n != ".DS_Store"}
    on_disk.discard("README.md")          # the manifest cannot hash itself

    assert not (on_disk - listed), f"packaged but unhashed: {sorted(on_disk - listed)}"
    assert not (listed - on_disk), f"manifest lists a file that is not packaged: {sorted(listed - on_disk)}"


def test_submission_manifest_sizes_and_hashes_are_current():
    """Every manifest row must still describe the file on disk.

    Separate from the name check above so the two failure modes stay distinguishable: a missing
    or extra filename and a stale digest have different causes and different fixes, and folding
    them together would let the first assertion mask the second. Presence in either direction is
    that test's business; this one only checks rows whose file is actually there.

    Mirrors tools/write_submission_readme.py: size from os.path.getsize rendered with `{:,}`,
    digest the first 16 hex characters of the file's SHA-256 (sha256() there, sliced at the row).
    A failure here means the manifest was written by an earlier build than the files beside it,
    which is the state a hand-run tools/package_supplementary_tables.py leaves behind.
    """
    sub = os.path.join(ROOT, "final_outputs")
    readme = os.path.join(sub, "README.md")
    if not os.path.isfile(readme):
        pytest.skip("final_outputs/ not built; run `make final-outputs`")

    with open(readme, encoding="utf-8") as fh:
        rows = re.findall(r"^\| `([^`]+)` \| ([\d,]+) \| `([0-9a-f]{16})` \|$", fh.read(), re.M)

    stale = []
    for name, size, digest in rows:
        path = os.path.join(sub, name)
        if not os.path.isfile(path):
            continue
        with open(path, "rb") as fh:
            actual = hashlib.sha256(fh.read()).hexdigest()[:16]
        on_disk_size = os.path.getsize(path)
        if actual != digest or on_disk_size != int(size.replace(",", "")):
            stale.append(f"{name}: manifest {int(size.replace(',', '')):,}B/{digest}, "
                         f"disk {on_disk_size:,}B/{actual}")

    assert not stale, (
        f"{len(stale)} of {len(rows)} manifest rows no longer describe the packaged file.\n  "
        + "\n  ".join(stale[:5])
        + (f"\n  ... and {len(stale) - 5} more" if len(stale) > 5 else "")
        + "\nThe manifest was written by an earlier build than the files beside it. "
          "Rerun `make final-outputs`, which rebuilds the tree and the manifest together.")


# --- GO over-representation correction family -------------------------------------

GO_PRIMARY = "figure_data/discordance_go_enrichment.csv"
GO_FULL = "figure_data/discordance_go_enrichment_full.csv"


def _go_rows(rel):
    path = os.path.join(ROOT, rel)
    # The full table is gitignored (written only by `make discordance-network`), so a fresh clone lacks it.
    if rel == GO_FULL and not os.path.isfile(path):
        pytest.skip("%s not built; run `make discordance-network`" % rel)
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("rel", [GO_PRIMARY, GO_FULL])
def test_go_family_size_is_recorded(rel):
    rows = _go_rows(rel)
    assert rows, "%s is empty" % rel
    assert "n_family" in rows[0], (
        "%s has no n_family column; the BH family size must be recorded for the fdr "
        "column to be interpretable" % rel)
    # BH runs once per (direction, background) block; the full table carries both backgrounds.
    per_block = {}
    for r in rows:
        block = (r["direction"], r.get("background", "coverage_restricted"))
        per_block.setdefault(block, set()).add(int(r["n_family"]))
    for block, sizes in per_block.items():
        assert len(sizes) == 1, (
            "%s: %s reports more than one family size %s" % (rel, "/".join(block), sorted(sizes)))


def test_go_bh_family_is_not_the_displayed_set():
    """BH must run over the pre-specified term family, not the rows that survive display.

    Selecting terms on k_fg and then correcting over only those terms makes the
    correction family a function of the test statistic. This is the defect the test exists
    to prevent, so it asserts the family is strictly larger than what is shown.
    """
    rows = _go_rows(GO_PRIMARY)
    for direction in sorted({r["direction"] for r in rows}):
        shown = [r for r in rows if r["direction"] == direction]
        family = int(shown[0]["n_family"])
        assert family > len(shown), (
            "%s: BH family is %d with %d rows displayed; the correction family must not "
            "equal the displayed set" % (direction, family, len(shown)))
        assert family >= 100, (
            "%s: a BH family of %d is too small for a term-size filtered GO family; check "
            "that zero-hit terms were retained" % (direction, family))


def test_go_full_family_retains_zero_hit_terms():
    """Zero-hit terms carry p = 1 and must stay in the family.

    Fisher's exact test is discrete; dropping the terms with no foreground gene removes
    the null distribution's atom at p = 1, leaving the retained nulls stochastically
    smaller than uniform, which makes BH anti-conservative.
    """
    full = _go_rows(GO_FULL)
    assert any(int(r["k_fg"]) == 0 for r in full), (
        "%s holds no zero-hit term, so the family was filtered on the foreground count"
        % GO_FULL)


def test_go_primary_sheet_holds_only_displayed_rows():
    for r in _go_rows(GO_PRIMARY):
        assert int(r["k_fg"]) >= 3, (
            "%s carries a row with k_fg=%s; the primary sheet is the display subset"
            % (GO_PRIMARY, r["k_fg"]))
