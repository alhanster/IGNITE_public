#!/usr/bin/env bash
# Per-step runner and logger.
#
# Usage: tools/run_step.sh <py|R> <script> [args...]
# Reads PY and RSCRIPT from the environment, exported by the Makefile. R_LIBS_USER is
# resolved here as well, so a direct call renders against the same pinned R stack.
# Writes one log per step to logs/NN_<script>.log and appends a record to logs/_summary.tsv.
# On non-zero exit, prints the failing target and the last 40 log lines to stderr, then propagates the original exit code.

set -uo pipefail

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <py|R> <script> [args...]" >&2
  exit 2
fi

KIND="$1"
TARGET="$2"
shift 2

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Resolve the pinned R stack here, not only in the Makefile. `make figures` exports
# R_LIBS_USER, but a direct `tools/run_step.sh R ...` call would otherwise fall back to the
# system library -- a different ragg renders byte-different PNGs from identical pixels, which
# silently splits tools/final_plots.sha256 across two encoders. `:=` mirrors the Makefile's `?=`:
# an R_LIBS_USER already in the environment wins.
if [ -d "$REPO_ROOT/.rlib" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-$REPO_ROOT/.rlib}"
fi

LOGDIR="$REPO_ROOT/logs"
SUMMARY="$LOGDIR/_summary.tsv"
mkdir -p "$LOGDIR"
[ -f "$SUMMARY" ] || printf 'step\tkind\ttarget\texit\tseconds\n' > "$SUMMARY"

# The header line counts as step 0, so the first real step is numbered 1.
STEP="$(wc -l < "$SUMMARY" | tr -d '[:space:]')"

TAG="$(printf '%02d' "$STEP")_$(basename "$TARGET" | tr -c 'A-Za-z0-9._-' '_' | sed 's/_*$//')"
LOG="$LOGDIR/$TAG.log"

t0=$SECONDS
rc=0
echo "==> [$(printf '%02d' "$STEP")] $TARGET"

case "$KIND" in
  py) "${PY:-python3}"      "$TARGET" "$@" > "$LOG" 2>&1 || rc=$? ;;
  R)  "${RSCRIPT:-Rscript}" "$TARGET" "$@" > "$LOG" 2>&1 || rc=$? ;;
  *)  echo "unknown step kind: $KIND (expected py or R)" >&2; exit 2 ;;
esac

printf '%s\t%s\t%s\t%s\t%s\n' \
  "$STEP" "$KIND" "$TARGET" "$rc" "$(( SECONDS - t0 ))" >> "$SUMMARY"

if [ "$rc" -ne 0 ]; then
  echo "FAILED (exit $rc): $TARGET" >&2
  echo "--- tail of $LOG ---" >&2
  tail -40 "$LOG" >&2
  exit "$rc"
fi
