#!/usr/bin/env bash
# Prints the best Python interpreter for this repo, or nothing if none qualifies.
# See REPRODUCIBILITY.md, Python (stage 1, `make tables`).
#
# Usage: tools/detect_python.sh   # prints e.g. "python3.12", or nothing
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REQ="$ROOT/requirements.txt"
[ -f "$REQ" ] || exit 0

# Unanchored on the leading '#' for Makefile copy-paste; matches only the line with both keyword and '>='.
MIN="$(sed -n 's/^.*python_requires: *>= *\([0-9][0-9.]*\).*/\1/p' "$REQ" | head -1)"
[ -n "$MIN" ] || exit 0

PIN="$(sed -n 's/^xgboost *== *\([0-9][^ ]*\).*/\1/p' "$REQ" | head -1)"

# .venv takes precedence over any python3 resolved from PATH, since it is provisioned by make setup.
CANDIDATES="$ROOT/.venv/bin/python3 python3 python3.12 python3.13 python3.14"

meets_floor() {
  "$1" -c 'import sys
m = tuple(int(x) for x in sys.argv[1].split("."))
sys.exit(0 if sys.version_info[:len(m)] >= m else 1)' "$MIN" 2>/dev/null
}

has_pin() {
  [ -n "$PIN" ] || return 1
  "$1" -c 'import sys
try:
    from importlib.metadata import version
    sys.exit(0 if version("xgboost") == sys.argv[1] else 1)
except Exception:
    sys.exit(1)' "$PIN" 2>/dev/null
}

# Pass 1: version floor met and exact pin installed, i.e. a ready-to-run environment.
for c in $CANDIDATES; do
  command -v "$c" >/dev/null 2>&1 || continue
  if meets_floor "$c" && has_pin "$c"; then echo "$c"; exit 0; fi
done

# Pass 2: floor met only, not runnable yet, but the guard's package error and fix apply to this interpreter.
for c in $CANDIDATES; do
  command -v "$c" >/dev/null 2>&1 || continue
  if meets_floor "$c"; then echo "$c"; exit 0; fi
done

exit 0
