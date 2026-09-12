#!/usr/bin/env python3
"""Fails at import time when installed package versions differ from the pins used to build committed artifacts. Called from src/analysis/pu_target_model.py via check_pins(). See REPRODUCIBILITY.md, Python (stage 1, `make tables`), for the numerical rationale and the IGNITE_ALLOW_VERSION_MISMATCH escape hatch."""
import os
import re
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REQUIREMENTS = os.path.join(ROOT, "requirements.txt")

ENV_OVERRIDE = "IGNITE_ALLOW_VERSION_MISMATCH"

# Packages whose exact version changes committed gene-level numbers. Only these are
# hard-enforced. Every other exact pin in requirements.txt is checked too, but a mismatch
# there is reported as a warning rather than an error (see ADVISORY below).
CRITICAL = ("xgboost",)

# Exact pins that are checked and reported but do not abort the run. A mismatch here has not
# been shown to move committed numbers, but it does mean the environment differs from the one
# the artifacts were built under, which is worth knowing before interpreting a diff.
ADVISORY = ("numpy", "pandas", "scipy", "scikit-learn", "pyarrow", "matplotlib")


class VersionPinError(RuntimeError):
    """Installed package version does not match the committed pin."""


class PythonVersionError(RuntimeError):
    """Interpreter is older than the pinned floor; kept separate from VersionPinError since the two need different fixes."""


def _parse_python_requires(path=REQUIREMENTS):
    """Returns the minimum Python version as a tuple, or None when requirements.txt has no python_requires directive, which skips the interpreter-floor check."""
    if not os.path.isfile(path):
        return None
    for raw in open(path):
        m = re.match(r"^#\s*python_requires:\s*>=\s*([0-9]+(?:\.[0-9]+)*)", raw.strip())
        if m:
            return tuple(int(x) for x in m.group(1).split("."))
    return None


def check_python(requirements=REQUIREMENTS):
    """Raises PythonVersionError if the interpreter is below the declared floor; runs before the package check, since a too-old interpreter could not have the pinned xgboost installed and the package-mismatch error would point at the wrong cause."""
    want = _parse_python_requires(requirements)
    if want is None:
        return None
    have = sys.version_info[:len(want)]
    if have >= want:
        return None

    w = ".".join(str(x) for x in want)
    h = ".".join(str(x) for x in have)
    msg = (
        "\n"
        + "=" * 72 + "\n"
        "WRONG PYTHON -- this is an interpreter problem, not a package problem\n"
        + "=" * 72 + "\n"
        "    required   Python >= %s\n"
        "    running    Python %s\n"
        "    from       %s\n\n"
        "The pinned xgboost publishes wheels only for Python >= %s, so it cannot\n"
        "be installed here at all. If you try, pip will report the newest version\n"
        "AVAILABLE TO THIS INTERPRETER rather than the newest that exists -- which\n"
        "reads as though the pin names a release that was never published. It was.\n\n"
        "A `python3` that is a conda base environment is the usual way to land here.\n\n"
        "Fix: point the build at a newer interpreter, e.g.\n"
        "    make tables PY=python3.12\n"
        "or create the environment from scratch:\n"
        "    python3.12 -m venv .venv && . .venv/bin/activate\n"
        "    pip install -r requirements.txt\n\n"
        "The Makefile normally auto-detects a suitable interpreter; reaching this\n"
        "message means no candidate on PATH satisfies the floor.\n\n"
        "To proceed anyway (do NOT commit the outputs): set %s=1\n"
        % (w, h, sys.executable, w, ENV_OVERRIDE)
        + "=" * 72
    )

    if os.environ.get(ENV_OVERRIDE, "").strip() not in ("", "0", "false", "False"):
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        print("[_version_guard] proceeding on Python %s (%s set)" % (h, ENV_OVERRIDE),
              file=sys.stderr)
        return (want, have)

    raise PythonVersionError(msg)


def _parse_exact_pins(path=REQUIREMENTS):
    """Return {package: version} for `pkg==X.Y.Z` lines. Comments stripped."""
    pins = {}
    if not os.path.isfile(path):
        return pins
    for raw in open(path):
        line = raw.split("#", 1)[0].strip()
        m = re.match(r"^([A-Za-z0-9_.\-]+)\s*==\s*([0-9][^\s,;]*)$", line)
        if m:
            pins[m.group(1).lower().replace("_", "-")] = m.group(2)
    return pins


def _installed(pkg):
    """Installed version string, or None if the package is absent."""
    try:
        from importlib.metadata import version, PackageNotFoundError
    except ImportError:                                    # py<3.8
        return None
    try:
        return version(pkg)
    except PackageNotFoundError:
        return None


def check_pins(critical=CRITICAL, requirements=REQUIREMENTS):
    """Raises VersionPinError on package pin mismatches; checks the interpreter first, since an incompatible interpreter makes the pinned xgboost uninstallable and would misdirect any mismatch found below it."""
    check_python(requirements)

    pins = _parse_exact_pins(requirements)
    mismatches = {}
    for pkg in critical:
        want = pins.get(pkg.lower())
        if want is None:
            continue                                       # not exactly pinned
        have = _installed(pkg)
        if have is None:
            # Missing critical packages are reported as failures rather than skipped, so `make check-versions` cannot report all pins clean on an interpreter with no xgboost installed.
            mismatches[pkg] = (want, "NOT INSTALLED")
            continue
        if have != want:
            mismatches[pkg] = (want, have)

    advisory = {}
    for pkg in ADVISORY:
        want = pins.get(pkg.lower())
        if want is None:
            continue
        have = _installed(pkg)
        if have is not None and have != want:
            advisory[pkg] = (want, have)
    if advisory:
        print("[_version_guard] environment differs from the build manifest in %d advisory "
              "pin(s); this does not abort the run:" % len(advisory), file=sys.stderr)
        for p, (w, h) in sorted(advisory.items()):
            print("    %-14s pinned %-10s installed %s" % (p, w, h), file=sys.stderr)

    if not mismatches:
        return mismatches

    detail = "\n".join(
        "    %-12s pinned %-10s installed %s" % (p, w, h)
        for p, (w, h) in sorted(mismatches.items()))
    msg = (
        "\n"
        + "=" * 72 + "\n"
        "VERSION PIN MISMATCH -- committed gene-level numbers will NOT reproduce\n"
        + "=" * 72 + "\n"
        + detail + "\n\n"
        "Committed artifacts (full_model_pu_scores.csv, panelD_recovery.csv,\n"
        "vignette_meta.json and the figures built from them) were produced under\n"
        "the pinned versions. A mismatched xgboost changes scores by ~2e-03,\n"
        "which is larger than the score gaps at the top of the ranking\n"
        "(0.00007-0.0027) -- so STAT4's rank and the panel-e recovery counts will\n"
        "differ. The run will otherwise look completely normal.\n\n"
        "Fix: run in the pinned environment, e.g.\n"
        "    make tables PY=/path/to/pinned/bin/python3\n"
        "or  pip install -r requirements.txt\n\n"
        "To proceed anyway (version-sensitivity work only -- do NOT commit the\n"
        "outputs): set %s=1\n" % ENV_OVERRIDE
        + "=" * 72
    )

    if os.environ.get(ENV_OVERRIDE, "").strip() not in ("", "0", "false", "False"):
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        print("[_version_guard] proceeding despite mismatch (%s set)" % ENV_OVERRIDE,
              file=sys.stderr)
        return mismatches

    raise VersionPinError(msg)


def audit_coverage(root=ROOT):
    """See REPRODUCIBILITY.md, `make check-versions`."""
    import glob
    unguarded = []
    pats = [os.path.join(root, "src", "**", "*.py"),
            os.path.join(root, "plots", "**", "*.py")]
    for pat in pats:
        for f in glob.glob(pat, recursive=True):
            if "__pycache__" in f or os.path.basename(f) == "_version_guard.py":
                continue
            src = open(f, encoding="utf-8", errors="replace").read()
            if "XGBClassifier(" not in src:
                continue
            # Counts a direct call to check_pins() or an inherited call via importing a module that already calls it.
            if "check_pins()" in src:
                continue
            if re.search(r"import\s+(attribution_decomposition|pu_target_model)"
                         r"|spec_from_file_location\(\s*[\"\']pum[\"\']", src):
                continue
            unguarded.append(os.path.relpath(f, root))
    return sorted(unguarded)


if __name__ == "__main__":
    # Prints the diagnosis without a traceback when run as a script; import-time callers still receive the exception and its traceback.
    try:
        bad = check_pins()
    except (PythonVersionError, VersionPinError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
    if bad:
        print("mismatch")
        sys.exit(1)
    gaps = audit_coverage()
    if gaps:
        print("version guard: pins OK, but these trainers are UNGUARDED:")
        for g in gaps:
            print("   ", g)
        print("\n  Add the guard snippet after their xgboost import "
              "(see any guarded script).")
        sys.exit(1)
    print("version guard: OK -- all critical pins match; "
          "all XGBoost trainers guarded")
