#!/usr/bin/env bash
# Regenerates final_plots.sha256 from the final_plots/ tree.
# Usage: tools/write_baseline.sh <py-interpreter> <rscript-interpreter>
# See REPRODUCIBILITY.md, Baseline generation, for prerequisites and hashing method.

set -euo pipefail

PY_BIN="${1:-}"
RS_BIN="${2:-}"
if [ -z "$PY_BIN" ] || [ -z "$RS_BIN" ]; then
  echo "usage: $0 <py-interpreter> <rscript-interpreter>" >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Resolves R packages from .rlib/, matching the environment that built the figures.
# See REPRODUCIBILITY.md, Environment layout and overrides.
if [ -d "$REPO_ROOT/.rlib" ]; then
  export R_LIBS_USER="$REPO_ROOT/.rlib"
fi

# Confirms the working tree matches HEAD before the baseline is written.
# See REPRODUCIBILITY.md, What a green gate does not cover.
DIRTY_PATHS="$(git status --porcelain -- \
                 src Makefile requirements.txt R-requirements.txt tools data figure_data \
               2>/dev/null || true)"
if [ -n "$DIRTY_PATHS" ]; then
  TREE_OK=no
else
  TREE_OK=yes
fi

XGB="$("$PY_BIN" -c 'import xgboost;print(xgboost.__version__)')"
RVER="$("$RS_BIN" -e 'cat(as.character(getRversion()))')"
GGP="$("$RS_BIN" -e 'cat(as.character(packageVersion("ggplot2")))')"
STAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Delegates toolchain-version checks to the existing checkers, which exit non-zero on mismatch.
# See REPRODUCIBILITY.md, What a green gate does not cover.
PINS_OK=yes
"$PY_BIN" src/analysis/common/_version_guard.py >/dev/null 2>&1 || PINS_OK=no
"$RS_BIN" tools/check_r_versions.R       >/dev/null 2>&1 || PINS_OK=no

N_FILES="$(find final_plots -type f ! -name .gitkeep ! -name .DS_Store -print0 | tr -dc '\0' | wc -c | tr -d ' ')"

# File classes come from git status, not a hardcoded list: regenerated=gitignored output, docs=tracked .md, frozen=other tracked.
N_REGEN=0; N_FROZEN=0; N_DOCS=0
while IFS= read -r -d '' f; do
  if git check-ignore -q "$f"; then
    N_REGEN=$((N_REGEN + 1))
  elif [ "${f##*.}" = "md" ]; then
    N_DOCS=$((N_DOCS + 1))
  else
    N_FROZEN=$((N_FROZEN + 1))
  fi
done < <(find final_plots -type f ! -name .gitkeep ! -name .DS_Store -print0)

OUT=final_plots.sha256
{
  echo "# IGNITE final_plots output baseline"
  echo "#"
  echo "# date:     $STAMP"
  echo "# xgboost:  $XGB   (pinned in requirements.txt; enforced by src/analysis/common/_version_guard.py)"
  echo "# R:        $RVER with ggplot2 $GGP   (pinned in R-requirements.txt; checked by make check-versions)"
  echo "#"
  if [ "$PINS_OK" = no ]; then
    echo "# ############################################################################"
    echo "# ## PROVISIONAL BASELINE -- THE TOOLCHAIN ABOVE DOES NOT MATCH THE PINS.   ##"
    echo "# ##                                                                        ##"
    echo "# ## The two version lines above describe the machine that ran this script, ##"
    echo "# ## and it is NOT the pinned environment. The hashes below may still be     ##"
    echo "# ## correct -- if the tree was rendered elsewhere under the pins and only   ##"
    echo "# ## the baseline was regenerated here -- but nothing in this file proves    ##"
    echo "# ## that. Do not treat a passing \`make verify\` against this baseline as     ##"
    echo "# ## evidence of reproducibility, and do not submit against it.              ##"
    echo "# ##                                                                        ##"
    echo "# ## Regenerate from a clean clone under the pinned interpreters before use. ##"
    echo "# ## Run \`make check-versions\` to see which packages differ.                 ##"
    echo "# ############################################################################"
    echo "#"
  fi
  if [ "$TREE_OK" = no ]; then
    echo "# ############################################################################"
    echo "# ## PROVISIONAL BASELINE -- HASHED FROM A DIRTY WORKING TREE.              ##"
    echo "# ##                                                                        ##"
    echo "# ## The commit line above is suffixed '-dirty' because build-relevant paths ##"
    echo "# ## (src/, Makefile, requirements, tools/, data/, figure_data/) held         ##"
    echo "# ## uncommitted changes when this file was written. The hashes describe the  ##"
    echo "# ## WORKING TREE, so the commit named above CANNOT reproduce them: checking  ##"
    echo "# ## it out and rebuilding yields different bytes for any figure the          ##"
    echo "# ## uncommitted work touches.                                                ##"
    echo "# ##                                                                        ##"
    echo "# ## Commit the work, rebuild from a clean clone, and regenerate before       ##"
    echo "# ## submitting or citing this baseline. Uncommitted paths at write time:     ##"
    printf '%s\n' "$DIRTY_PATHS" | sed 's/^/# ##   /'
    echo "# ############################################################################"
    echo "#"
  fi
  echo "# produced: git clone <IGNITE> /tmp/verify && cd /tmp/verify && \\"
  echo "#             make tables PY=<pinned-python> RSCRIPT=<pinned-rscript> && \\"
  echo "#             make figures PY=<pinned-python> RSCRIPT=<pinned-rscript> && \\"
  echo "#             tools/write_baseline.sh <pinned-python> <pinned-rscript>"
  echo "#"
  echo "# IMPORTANT: the regenerated files hashed below are BUILD OUTPUTS. They are"
  echo "# gitignored and are NOT in the commit -- you cannot 'git checkout' them. Rebuild"
  echo "# with \`make figures\` (fast, no model fitting) and check with \`make verify\`."
  echo "#"
  echo "# THE R VERSION ABOVE IS LOAD-BEARING. Rendering libraries change pixel output"
  echo "# between releases and nothing guards R at render time, so a mismatched R stack"
  echo "# produces figure mismatches that look like code regressions. Run"
  echo "# \`make check-versions\` first; it checks both stacks."
  echo "#"
  echo "# WHAT THIS FILE DOES AND DOES NOT PROVE. Under the pinned R stack recorded above,"
  echo "# the PNG hashes DO reproduce: 10 of 10 were byte-identical across two independent"
  echo "# clones rendered through \`make figures\`. So a PNG mismatch on this machine is a real"
  echo "# signal, usually a render that did not resolve .rlib and used a different ragg."
  echo "# What is not pinnable is the native rasterization stack -- ragg goes through"
  echo "# freetype/libpng and the installed fonts (see R-requirements.txt) -- so a DIFFERENT"
  echo "# MACHINE can reproduce every NUMBER exactly and still fail these hashes."
  echo "#"
  echo "# So \`make verify-figures\` reports rendered-PNG differences without failing, and the"
  echo "# portable claim lives in \`make verify-tables\`: figure_data/ byte-identical to the"
  echo "# commit after a rebuild. That is the check to trust for reproduction; this one is a"
  echo "# same-machine regression check plus a tamper check on the artifacts no script can"
  echo "# rebuild."
  echo "#"
  echo "# $N_FILES files: $N_REGEN regenerated outputs (the reproducibility-relevant subset),"
  echo "# $N_FROZEN frozen artifacts with no producing script, and $N_DOCS figure-provenance docs."
  echo "# Those three counts are COMPUTED, not typed."
  echo "#"
  echo "# Was 38 before the two-stage restructure. What left, and why:"
  echo "#   5  internal working documents, moved out of the repository entirely"
  echo "#   1  causal_functional_genomics.png, dropped with its renderer as a strict subset"
  echo "#      of the _with_thde figure"
  echo "#   1  the pipeline schematic, DELETED -- stale depth=3 artwork, no producing script"
  echo "#"
  find final_plots -type f ! -name .gitkeep ! -name .DS_Store -print0 \
    | sort -z \
    | xargs -0 shasum -a 256
} > "$OUT"

echo "wrote $OUT"
echo "  $(grep -cE '^[0-9a-f]{64}  ' "$OUT") hashed files"
echo "  xgboost $XGB / R $RVER + ggplot2 $GGP"
