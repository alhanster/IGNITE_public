"""Whether this run is on the reference platform, where committed tables reproduce exactly.

The committed figure_data/ tables were produced on Apple Silicon macOS. Elsewhere XGBoost and BLAS
floating point drift, so checks that pin a value to the committed run warn instead of failing.
See REPRODUCIBILITY.md.
"""
import platform


def on_reference_platform():
    return platform.system() == "Darwin" and platform.machine() == "arm64"


def check(ok, msg):
    """Fails on the reference platform; elsewhere prints msg as a warning and continues."""
    if ok:
        return
    if on_reference_platform():
        raise AssertionError(msg)
    print(f"  WARNING (not the reference platform): {msg}")
