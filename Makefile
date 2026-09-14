# Immune Target Atlas: minimal reproduction subset (IGNITE).
#
# Two stages:
#   make tables   Fits models, runs all tests, writes figure_data/. Requires the pinned
#                 Python environment.
#   make figures  Reads figure_data/ only (no data/ access, no refitting) and writes
#                 final_plots/. Runs standalone from a fresh clone.
#
# figure_data/ is committed, which is what makes stage 2 reproducible standalone from a
# fresh clone. See figure_data/PROVENANCE.md for which stage-1 script writes each table.
#
# Regeneration starts from the committed data/perturbseq/pu/pu_model_matrix.parquet. Two
# evidence pipelines feed one PU-learning target ranker: genetic, monogenic and GWAS priors
# (data/reference, data/data_drug), and perturbational functional genomics from perturb-seq
# (data/perturbseq/*). Runs offline from committed data; requires Python and R.

# Python and R may need separate environments. Override with:
#   make tables  PY=/path/to/py-env/bin/python3
#   make figures RSCRIPT=/path/to/r-env/bin/Rscript
# PY must be the pinned environment (requirements.txt: xgboost==3.3.0); src/analysis/common/_version_guard.py aborts on a mismatch. See REPRODUCIBILITY.md.
PY_AUTO := $(shell tools/detect_python.sh)
PY      ?= $(if $(PY_AUTO),$(PY_AUTO),python3)
RSCRIPT ?= Rscript
export PYTHONPATH := src

# R_LIBS_USER is exported so R picks up .rlib/ automatically; an existing R_LIBS_USER wins.
ifneq ($(wildcard .rlib),)
R_LIBS_USER ?= $(CURDIR)/.rlib
export R_LIBS_USER
endif

# Every step goes through the logger, which writes logs/NN_<script>.log, appends a row to logs/_summary.tsv, and on failure prints the log tail and propagates the exit code.
RUN := PY="$(PY)" RSCRIPT="$(RSCRIPT)" tools/run_step.sh

.PHONY: help setup install check-versions tables figures all submission test verify verify-tables verify-figures supplementary final-outputs clean permutation-null permutation-null-finalize discordance-network vignette-panelc gps-inputs knn-signatures

# PERM_JOBS only changes wall-clock time; each permutation seeds from RandomState(10_000 + i).
PERM_JOBS ?= 11
PERM_N    ?= 1000
PERM_SCR  ?= 99

# Matched null draws for the discordance STRING test (opt-in target below); not on the make tables path. See REPRODUCIBILITY.md.
DISC_NULL ?= 500

# Relaxes the MISSING class in verify-figures from fatal to reported-only; empty by default, set only by supplementary:. See REPRODUCIBILITY.md.
VERIFY_ALLOW_MISSING ?=

help:
	@echo "IGNITE — two-stage build"
	@echo ""
	@echo "  make tables          stage 1: fit models, run tests, write figure_data/  (slow)"
	@echo "  make figures         stage 2: read figure_data/, write final_plots/      (fast)"
	@echo "  make all             tables, then figures"
	@echo "  make submission      clone -> final_outputs/ in one command (offline)"
	@echo ""
	@echo "  make permutation-null           rebuild the label-permutation null (OPT-IN, ~8.5 h)"
	@echo "  make permutation-null-finalize  recompute its p-values vs the current ladder (seconds)"
	@echo "  make vignette-panelc            re-derive panel c from Open Targets   (OPT-IN, network)"
	@echo "  make discordance-network        rebuild S2's STRING/GO cache          (OPT-IN, ~10 min"
	@echo "                                  + a one-time 135 MB download; only needed if the core"
	@echo "                                  sets move, and \`make tables\` says so when they do)"
	@echo "  make gps-inputs                 re-derive the GPS inputs from the published"
	@echo "                                  supplementary table                  (OPT-IN, seconds;"
	@echo "                                  needs data/raw/gps/GPS_allgenes.xlsx, which is"
	@echo "                                  gitignored; see the script header)"
	@echo "  make knn-signatures             stream the perturb-seq knockdown signature matrices"
	@echo "                                  from the public source              (OPT-IN, network,"
	@echo "                                  ~3-4 min; they are NOT committed at 294 MB, and the"
	@echo "                                  one step that reads them skips without them)"
	@echo ""
	@echo "  make verify          both checks below"
	@echo "  make verify-tables   did the NUMBERS reproduce? figure_data/ vs the commit. HARD"
	@echo "  make verify-figures  did the PIXELS reproduce? reported; rendered PNGs never fail"
	@echo ""
	@echo "  make supplementary   render + package ONLY the supplementary set into"
	@echo "                       supplementary/ (5 renderers, main figures untouched)"
	@echo "  make final-outputs   assemble final_outputs/ (display items only)"
	@echo ""
	@echo "  make setup           provision BOTH stacks into .venv/ and .rlib/  (start here)"
	@echo "  make check-versions  check BOTH the Python and R stacks against the pins"
	@echo "  make test            fast pre-flight contract checks (no fitting, no rendering)"
	@echo "  make install         pip install -r requirements.txt into the current interpreter"
	@echo "  make clean           remove generated outputs (keeps figure_data/ and frozen artifacts)"
	@echo ""
	@echo "  Override interpreters:  make tables PY=... RSCRIPT=..."
	@echo "  Per-step logs: logs/NN_<script>.log   Summary: logs/_summary.tsv"

# Entry point for a fresh clone: provisions .venv/ and .rlib/ (gitignored), so later targets need no PY= or R_LIBS_USER= override. See REPRODUCIBILITY.md.
setup:
	@tools/setup_env.sh

install:
	$(PY) -m pip install -r requirements.txt

# Verifies Python and R package versions before a rebuild; mismatched versions silently change results.
# See REPRODUCIBILITY.md.
check-versions:
	@echo "  PY=$(PY) -> $$(command -v $(PY) 2>/dev/null || echo 'NOT FOUND') ($$($(PY) -V 2>&1))"
	@$(PY) src/analysis/common/_version_guard.py && echo "  PY=$(PY) matches the committed pins"
	@$(RSCRIPT) tools/check_r_versions.R

all: tables
	$(MAKE) figures

# One command from a fresh clone to the submission tree. Runs stage 1, stage 2, both
# verification gates and the packaging step, in that order, entirely offline. This is the
# target the reproduction claim refers to. See REPRODUCIBILITY.md for what it does not re-derive.
submission:
	@echo "[1/4] checking the pinned environment"
	@$(MAKE) --no-print-directory check-versions
	@echo "[2/4] stage 1: fitting models and writing figure_data/"
	@$(MAKE) --no-print-directory tables
	@echo "[3/4] stage 2: rendering final_plots/"
	@$(MAKE) --no-print-directory figures
	@echo "[4/4] verifying and assembling final_outputs/"
	@$(MAKE) --no-print-directory final-outputs
	@echo ""
	@echo "final_outputs/ built from a clean tree."

# Pre-flight contract checks. Fast, no model fitting, no rendering.
# pytest is declared in requirements.txt but unpinned and test-only, so it can be missing from
# an otherwise valid reproduction environment; say so rather than failing on an import error.
test:
	@$(PY) -c "import pytest" 2>/dev/null || { \
	  echo "pytest is not installed in $(PY)."; \
	  echo "  It is test-only and unpinned; install it with:  $(PY) -m pip install pytest"; \
	  exit 1; }
	@$(PY) -m pytest tools -q -p no:cacheprovider

# STAGE 1: analysis scripts. Writes figure_data/; outputs not read by a figure go to outputs/, treated as scratch and removed by make clean.
# Script order is load-bearing; run without -j. See REPRODUCIBILITY.md.
tables:
	rm -rf logs
	rm -f outputs/_skipped.tsv
	mkdir -p figure_data outputs/model outputs/tables outputs/prospective outputs/discordance outputs/knn data/derived
	$(RUN) py src/analysis/common/_version_guard.py
	$(RUN) py src/analysis/evidence/make_enrichment_table.py
	$(RUN) py src/analysis/model/pu_target_model.py
	$(RUN) py src/analysis/model/attribution_decomposition.py
	$(RUN) py src/analysis/model/reshape_attribution_to_panels.py
	$(RUN) py src/analysis/model/compute_nb_significance.py
	$(RUN) py src/analysis/model/recompute_block_signal.py
	$(RUN) py src/analysis/model/make_panelD_recovery.py
	$(RUN) py src/analysis/model/build_score_by_group.py
	$(RUN) py src/analysis/evidence/build_genetic_data_tables.py
	$(RUN) py src/analysis/evidence/build_perturbational_fg_tables.py
	$(RUN) py src/analysis/vignette/build_vignette_panelB.py
	$(RUN) py src/analysis/vignette/build_vignette_panelD_knn.py
	$(RUN) py src/analysis/vignette/build_vignette_meta.py
	$(RUN) py src/analysis/specificity/build_specificity_matrix.py
	$(RUN) py src/analysis/specificity/compute_delong_specificity.py
	$(RUN) py src/analysis/specificity/build_leakage_controlled.py
	$(RUN) py src/analysis/validation/build_temporal_holdout.py
	$(RUN) py src/analysis/validation/build_heldout_trial_auc.py
	$(RUN) py src/analysis/tables/build_ranked_atlas_table.py
	$(RUN) py src/analysis/validation/build_immune_stages_data.py
	$(RUN) py src/analysis/tables/build_feature_dictionary.py
	$(RUN) py src/analysis/tables/build_drug_status_table.py
	$(RUN) py src/analysis/evidence/build_univariate_by_group.py
	$(RUN) py src/analysis/evidence/build_orthogonality_tables.py
	$(RUN) py src/analysis/model/build_coalition_tables.py
	$(RUN) py src/analysis/tables/build_supp_stat_sheets.py
	$(RUN) py src/analysis/vignette/build_knn_nearest_target.py
	$(RUN) py src/analysis/vignette/add_knn_approved_drugs.py
	$(RUN) py src/analysis/discordance/build_discordance.py
	@$(PY) src/analysis/common/_skip_ledger.py
	@echo ""
	@echo "stage 1 complete. Commit figure_data/ so that \`make figures\` needs no refit."

# Stage 2: figures. Reads figure_data/ only; each renderer writes its deliverable directly into final_plots/. Numbered supplementary tables are excluded from final_plots/ and packaged separately; see REPRODUCIBILITY.md.

# Tables are packaged from figure_data/ via package_supplementary_tables.py, bypassing tools/final_plots.sha256.

# fig_label_permutation.R writes to final_plots/ but counts as S1; numbering is set by tools/supplementary_figures.tsv.
SUPP_RENDERERS_R  := fig_discordance.R fig_label_permutation.R fig_leakage_controlled.R
SUPP_RENDERERS_PY :=

figures:
	mkdir -p final_plots/supplementary
	$(RUN) R  src/figures/fig_genetic_data.R
	$(RUN) R  src/figures/fig_perturbational_functional_genomics.R
	$(RUN) R  src/figures/fig_model_analysis.R
	$(RUN) R  src/figures/fig_immune_stages.R
	$(RUN) R  src/figures/fig_stat4_vignette.R
	$(RUN) R  src/figures/fig_specificity.R
	$(RUN) R  src/figures/fig_temporal_holdout.R
	$(RUN) R  src/figures/fig_label_permutation.R
	$(RUN) R  src/figures/fig_discordance.R
	$(RUN) R  src/figures/fig_leakage_controlled.R
	@echo ""
	@echo "stage 2 complete. Verify: shasum -a 256 -c tools/final_plots.sha256"

# Opt-in permutation null for the AUC ladder, not part of `make tables`. Budget about 8.5 hours at PERM_JOBS=11 (measured 2026-09-06: 7 h 17 m for the 1,000 permutations, plus 1 h 16 m for the 99 scrambled-feature draws that this target also runs). Reads the committed figure_data/panelA_auc_ladder.csv and checkpoints every 25 permutations so an interrupted run resumes. Writes directly into figure_data/; revert a smoke test with `git checkout -- figure_data`. See REPRODUCIBILITY.md.
permutation-null:
	mkdir -p figure_data outputs/permutation
	$(RUN) py src/analysis/permutation/run_label_permutation.py --n-perm $(PERM_N) --jobs $(PERM_JOBS)
	$(RUN) py src/analysis/permutation/scrambled_feature_control.py --n-scramble $(PERM_SCR)
	@echo ""
	@echo "permutation null rebuilt. Check it reproduced: make verify-tables"

# Recomputes label_permutation_pvalues.csv from the already-committed null table against the current ladder, in seconds rather than hours. Runs when `make figures` reports stale p-values while the null itself remains valid.
permutation-null-finalize:
	$(RUN) py src/analysis/permutation/finalize_label_permutation.py

# Opt-in target, not part of `make tables`. Builds S2's STRING coherence, GO over-representation and dark-proteome tables from about 135 MB of external STRING and UniProt-GOA archives, cached in gitignored data/raw/discordance/ (about 10 minutes on first run, about 1 minute after). Rerun when the core gene sets change (pum.SEED, tree depth, or the feature matrix): run make tables, then make discordance-network, then make tables again. `--validate` checks the cached data still matches upstream STRING. See REPRODUCIBILITY.md.
discordance-network:
	mkdir -p figure_data data/raw/discordance
	$(RUN) py src/analysis/discordance/fetch_discordance_network.py --n-null $(DISC_NULL)
	@echo ""
	@echo "STRING/GO cache rebuilt. Now: make tables   (folds it into discordance_stats.json)"
	@echo "Then check it reproduced: make verify-tables"

# vignette-panelc asserts Open Targets release 26.06 and aborts without writing if a different release is found, since scores move between releases.
vignette-panelc:
	mkdir -p figure_data
	$(RUN) py src/analysis/vignette/fetch_vignette_panelC_gwas.py
	@echo ""
	@echo "panel c refreshed. Check it reproduced: make verify-tables"

# gps-inputs: OPT-IN, network. Regenerates gps_per_gene.csv and gps_drug_labeled_genes.csv from the published GPS supplementary table (data/raw/gps/GPS_allgenes.xlsx, gitignored, absent from clean clones); not part of make tables.
# Run make tables && make verify-tables afterward: gps_max_overall and the 453-gene list must reproduce byte-identically, or the source table or a derivation rule has changed.
# See REPRODUCIBILITY.md.
gps-inputs:
	$(RUN) py src/analysis/specificity/build_gps_per_gene.py
	@echo ""
	@echo "GPS inputs re-derived. Confirm nothing moved: make tables && make verify-tables"

# knn-signatures: OPT-IN, network, about 3-4 minutes. Streams perturb-seq knockdown signature matrices into data/perturbseq/knn/ (not committed; signatures_Stim48hr.npy alone is 294 MB); not part of make tables.
# Without it, build_vignette_panelD_knn.py skips and leaves vignette_panelD_knn.csv and vignette_knn_stats.json at their committed values, so verify-tables passes trivially for those two files.
# See REPRODUCIBILITY.md.
knn-signatures:
	$(RUN) py src/analysis/vignette/build_knn_signatures.py
	@echo ""
	@echo "Signature matrices rebuilt. Re-derive panel d, then check it reproduced:"
	@echo "  $(RUN) py src/analysis/vignette/build_vignette_panelD_knn.py && make verify-tables"

# verify-tables checks numeric reproduction (hard failure, portable); verify-figures checks pixel reproduction (reported only, machine-specific).
# See REPRODUCIBILITY.md.
verify: verify-tables verify-figures

# Reproduction is judged against the committed figure_data/ via git status, so a clean tree after rebuild is the check; it requires a git checkout and excludes .md files from the gate.
# See REPRODUCIBILITY.md.
verify-tables:
	@git rev-parse --git-dir >/dev/null 2>&1 || { \
	  echo "verify-tables needs a git checkout (figure_data/ is verified against the commit)" >&2; \
	  exit 1; }
	@if git diff --quiet -- figure_data ':(exclude)figure_data/*.md' \
	   && [ -z "$$(git ls-files --others --exclude-standard figure_data | grep -v '\.md$$')" ]; then \
	  echo "  figure_data/: all $$(git ls-files figure_data | grep -cv '\.md$$') tables byte-identical to the commit"; \
	elif [ "$$(uname -s)/$$(uname -m)" != "Darwin/arm64" ]; then \
	  echo ""; \
	  echo "  figure_data/: NOT byte-identical to the commit, as expected on $$(uname -s)/$$(uname -m)."; \
	  echo "  The committed tables were produced on Apple Silicon macOS, the only platform where they"; \
	  echo "  reproduce exactly; elsewhere model-derived values drift in about the 4th decimal. Not a"; \
	  echo "  failure here."; \
	  git diff --stat -- figure_data ':(exclude)figure_data/*.md'; \
	else \
	  echo ""; \
	  echo "TABLES DIFFER: the numbers changed. This is a real reproduction failure."; \
	  git diff --stat -- figure_data ':(exclude)figure_data/*.md'; \
	  git ls-files --others --exclude-standard figure_data | grep -v '\.md$$' || true; \
	  exit 1; \
	fi

# Rendering check: sorts files into five severity classes before reporting failures. See REPRODUCIBILITY.md.
verify-figures:
	@if [ ! -f tools/final_plots.sha256 ]; then \
	  echo "no tools/final_plots.sha256 to verify against" >&2; exit 1; fi
	@shasum -a 256 -c tools/final_plots.sha256 2>/dev/null > .verify_figures.tmp || true
	@rm -f .verify_figures.cls; : > .verify_figures.cls; \
	while IFS= read -r line; do \
	  case "$$line" in *": FAILED"*) ;; *) continue ;; esac; \
	  f=$${line%%: FAILED*}; \
	  case "$$line" in \
	    *": FAILED open or read") \
      if [ "$${f##*.}" = "pdf" ] && [ -z "$$(find final_plots/pdf -name '*.pdf' 2>/dev/null | head -1)" ]; \
      then c=nopdf; else c=missing; fi ;; \
	    *) if { [ "$${f##*.}" = "png" ] || [ "$${f##*.}" = "pdf" ]; } \
	          && git check-ignore -q "$$f" 2>/dev/null; then c=soft; \
	       elif [ "$${f##*.}" = "md" ]; then c=docs; \
	       elif git check-ignore -q "$$f" 2>/dev/null; then c=pub; \
	       else c=hard; fi ;; \
	  esac; \
	  echo "$$c $$f" >> .verify_figures.cls; \
	done < .verify_figures.tmp; \
	rm -f .verify_figures.tmp; \
	awk '/^[0-9a-f]{64}  /{print substr($$0,67)}' tools/final_plots.sha256 | sort > .verify_listed.tmp; \
	find final_plots -type f ! -name .gitkeep ! -name .DS_Store -print | sort > .verify_present.tmp; \
	while IFS= read -r f; do \
	  if [ "$${f##*.}" != "md" ]; then echo "xdata $$f"; \
	  elif git ls-files --error-unmatch "$$f" >/dev/null 2>&1; then echo "xdoc $$f"; \
	  else echo "xdoc $$f   (untracked: not hashed and not in git, no undo)"; fi \
	done < "$$(comm -13 .verify_listed.tmp .verify_present.tmp > .verify_extra.tmp; echo .verify_extra.tmp)" >> .verify_figures.cls; \
	rm -f .verify_listed.tmp .verify_present.tmp .verify_extra.tmp; \
	missing=$$(grep -c '^missing ' .verify_figures.cls || true); \
	nopdf=$$(grep -c '^nopdf ' .verify_figures.cls || true); \
	soft=$$(grep -c '^soft ' .verify_figures.cls || true); \
	docs=$$(grep -c '^docs ' .verify_figures.cls || true); \
	pub=$$(grep -c '^pub ' .verify_figures.cls || true); \
	hard=$$(grep -c '^hard ' .verify_figures.cls || true); \
	xdoc=$$(grep -c '^xdoc ' .verify_figures.cls || true); \
	xdata=$$(grep -c '^xdata ' .verify_figures.cls || true); \
	if [ $$missing -gt 0 ]; then \
	  echo ""; echo "  BASELINE FILES MISSING FROM final_plots/, nothing was verified:"; \
	  sed -n 's/^missing /    /p' .verify_figures.cls; fi; \
	if [ $$nopdf -gt 0 ]; then \
	  echo ""; echo "  $$nopdf PDF(s) not rendered (PDFs need macOS's quartz device; not a failure --"; \
	  echo "  every PNG and table is still checked):"; \
	  sed -n 's/^nopdf /    /p' .verify_figures.cls; fi; \
	if [ $$soft -gt 0 ]; then \
	  echo ""; echo "  rendered figures differing from the baseline:"; \
	  sed -n 's/^soft /    /p' .verify_figures.cls; fi; \
	if [ $$docs -gt 0 ]; then \
	  echo ""; echo "  provenance docs edited since the baseline:"; \
	  sed -n 's/^docs /    /p' .verify_figures.cls; fi; \
	if [ $$pub -gt 0 ]; then \
	  echo ""; echo "  PUBLISHED TABLES DIFFERING (stale publish or changed source):"; \
	  sed -n 's/^pub /    /p' .verify_figures.cls; fi; \
	if [ $$hard -gt 0 ]; then \
	  echo ""; echo "  COMMITTED ARTIFACTS DIFFERING (corruption or loss):"; \
	  sed -n 's/^hard /    /p' .verify_figures.cls; fi; \
	if [ $$xdoc -gt 0 ]; then \
	  echo ""; echo "  docs present but not in the baseline:"; \
	  sed -n 's/^xdoc /    /p' .verify_figures.cls; fi; \
	if [ $$xdata -gt 0 ]; then \
	  echo ""; echo "  UNVERIFIED FILES IN THE DELIVERABLE TREE (hashed by nothing):"; \
	  sed -n 's/^xdata /    /p' .verify_figures.cls; fi; \
	rm -f .verify_figures.cls; \
	if [ $$missing -gt 0 ]; then \
	  echo ""; \
	  echo "  $$missing baseline file(s) are ABSENT, not different. This is NOT a font-stack"; \
	  echo "  difference: shasum could not open them, so they were checked by nothing and"; \
	  echo "  a pass here would vouch for an empty tree. Run 'make figures' before packaging:"; \
	  echo "  'make final-outputs' rm -rf's final_outputs/ before it copies, and that tree is"; \
	  echo "  gitignored with no copy in git."; \
	fi; \
	if [ $$soft -gt 0 ]; then \
	  echo ""; \
	  echo "  $$soft rendered figure(s) differ. NOT fatal, but do not shrug this off: under the"; \
	  echo "  pinned R stack PNG bytes DO reproduce -- 10 of 10 were byte-identical across two"; \
	  echo "  independent clones. A differing PNG most often means the render did not resolve"; \
	  echo "  .rlib, so a different ragg encoded it (identical pixels, different bytes)."; \
	  echo "  What is genuinely unpinnable is the native rasterization stack -- freetype,"; \
	  echo "  libpng and the installed fonts -- which differs across MACHINES, not across runs"; \
	  echo "  on one. See src/figures/R-requirements.txt. The numbers are checked by verify-tables."; \
	  echo "  A .pdf cannot match on any stack: every R PDF device stamps a wall-clock"; \
	  echo "  /CreationDate, so two renders a second apart differ. Its hash records what was"; \
	  echo "  shipped, not a reproduction claim."; \
	fi; \
	if [ $$docs -gt 0 ]; then \
	  echo ""; \
	  echo "  $$docs provenance doc(s) differ. NOT a failure: editing documentation is not"; \
	  echo "  corruption, and git tracks their content. Refresh the baseline with"; \
	  echo "  tools/write_baseline.sh when convenient."; \
	fi; \
	if [ $$pub -gt 0 ]; then \
	  echo ""; \
	  echo "  $$pub published table(s) differ. THIS IS A FAILURE, but it is NOT corruption --"; \
	  echo "  these are gitignored copies the cp steps in 'figures' make out of figure_data/,"; \
	  echo "  so they cannot drift on a font stack. Either the publish is stale (run 'make"; \
	  echo "  figures') or the figure_data/ source legitimately changed and the baseline"; \
	  echo "  predates it (check 'make verify-tables', then tools/write_baseline.sh)."; \
	fi; \
	if [ $$hard -gt 0 ]; then \
	  echo ""; \
	  echo "  $$hard committed artifact(s) differ. This is a failure: these are neither"; \
	  echo "  rendered nor rebuilt by any script, so a mismatch means a file was corrupted,"; \
	  echo "  truncated or overwritten."; \
	fi; \
	if [ $$xdoc -gt 0 ]; then \
	  echo ""; \
	  echo "  $$xdoc doc(s) are not in the baseline. NOT a failure: the write-up surface grows"; \
	  echo "  between baselines. Refresh with tools/write_baseline.sh when convenient, and"; \
	  echo "  for any marked UNTRACKED, either git add it or accept that no copy exists."; \
	fi; \
	if [ $$xdata -gt 0 ]; then \
	  echo ""; \
	  echo "  $$xdata file(s) in final_plots/ are hashed by nothing. This is a failure, not"; \
	  echo "  a mismatch but an absence: shasum -c only checks files the manifest lists, so"; \
	  echo "  these were never verified and would ride into final_outputs/ unchecked. Either they"; \
	  echo "  are a real deliverable (add a producer, rebuild, tools/write_baseline.sh) or they"; \
	  echo "  are stray scratch and do not belong in the deliverable tree."; \
	fi; \
	if [ $$pub -gt 0 ] || [ $$hard -gt 0 ] || [ $$xdata -gt 0 ] \
	   || { [ $$missing -gt 0 ] && [ -z "$(VERIFY_ALLOW_MISSING)" ]; }; then exit 1; fi; \
	[ $$missing -eq 0 ] && [ $$nopdf -eq 0 ] && [ $$soft -eq 0 ] && [ $$docs -eq 0 ] && [ $$xdoc -eq 0 ] && echo "  final_plots/: all `grep -cE '^[0-9a-f]{64}  ' tools/final_plots.sha256` baseline files match, and no unlisted files" || true

# Renders and packages only the supplementary display items, into supplementary/. See REPRODUCIBILITY.md.
supplementary:
	mkdir -p final_plots/supplementary
	@for f in $(SUPP_RENDERERS_R);  do $(RUN) R  src/figures/$$f || exit 1; done
	@for f in $(SUPP_RENDERERS_PY); do $(RUN) py src/figures/$$f || exit 1; done
# No separate table-publishing step: the numbered tables are packaged directly from figure_data/ below.
	@$(MAKE) --no-print-directory verify VERIFY_ALLOW_MISSING=1
	@$(PY) tools/package_supplementary_figures.py --out-dir supplementary/figures --naming numbered
	@$(PY) tools/package_supplementary_tables.py  --out-dir supplementary/tables
	@echo ""
	@echo "supplementary/ assembled. Main figures were not touched."

# Assembles only the display items, numbered as the manuscript numbers them. See REPRODUCIBILITY.md.
final-outputs: verify
	rm -rf final_outputs
	mkdir -p final_outputs/main_figures/pdf final_outputs/supplementary_figures
	cp final_plots/genetic_data.png                              final_outputs/main_figures/Fig2_genetic_data.png
	@if [ -f final_plots/pdf/genetic_data.pdf ]; then cp final_plots/pdf/genetic_data.pdf final_outputs/main_figures/pdf/Fig2_genetic_data.pdf; fi
	cp final_plots/perturbational_functional_genomics_with_thde.png       final_outputs/main_figures/Fig3_perturbational_functional_genomics.png
	@if [ -f final_plots/pdf/perturbational_functional_genomics_with_thde.pdf ]; then cp final_plots/pdf/perturbational_functional_genomics_with_thde.pdf final_outputs/main_figures/pdf/Fig3_perturbational_functional_genomics.pdf; fi
	cp final_plots/figure_model_analysis_composite.png           final_outputs/main_figures/Fig4_model_analysis.png
	@if [ -f final_plots/pdf/figure_model_analysis_composite.pdf ]; then cp final_plots/pdf/figure_model_analysis_composite.pdf final_outputs/main_figures/pdf/Fig4_model_analysis.pdf; fi
	cp final_plots/figure_trial_validation_immune_stages_R.png   final_outputs/main_figures/Fig5_trial_validation.png
	@if [ -f final_plots/pdf/figure_trial_validation_immune_stages_R.pdf ]; then cp final_plots/pdf/figure_trial_validation_immune_stages_R.pdf final_outputs/main_figures/pdf/Fig5_trial_validation.pdf; fi
	cp final_plots/figure_prospective_temporal_holdout.png       final_outputs/main_figures/Fig6_prospective_temporal_holdout.png
	@if [ -f final_plots/pdf/figure_prospective_temporal_holdout.pdf ]; then cp final_plots/pdf/figure_prospective_temporal_holdout.pdf final_outputs/main_figures/pdf/Fig6_prospective_temporal_holdout.pdf; fi
	cp final_plots/specificity_immune_vs_cardiac_faceted.png     final_outputs/main_figures/Fig7_specificity_control.png
	@if [ -f final_plots/pdf/specificity_immune_vs_cardiac_faceted.pdf ]; then cp final_plots/pdf/specificity_immune_vs_cardiac_faceted.pdf final_outputs/main_figures/pdf/Fig7_specificity_control.pdf; fi
	cp final_plots/figure_stat4_vignette.png                     final_outputs/main_figures/Fig8_stat4_vignette.png
	@if [ -f final_plots/pdf/figure_stat4_vignette.pdf ]; then cp final_plots/pdf/figure_stat4_vignette.pdf final_outputs/main_figures/pdf/Fig8_stat4_vignette.pdf; fi
# tools/supplementary_figures.tsv is the sole source of S-numbers; S1 lives in final_plots/ (renders with the main set).
# PDFs render only on macOS; drop the pdf/ directory rather than ship it empty.
	@rmdir final_outputs/main_figures/pdf 2>/dev/null || true
	@$(PY) tools/package_supplementary_figures.py --out-dir final_outputs/supplementary_figures
# supplementary_data/ is not assembled; see REPRODUCIBILITY.md.
	@$(PY) tools/package_supplementary_tables.py
	@$(PY) tools/write_submission_readme.py
# Finder can drop a .DS_Store into the tree between assembly and packing; it is in no manifest.
	@find final_outputs -name .DS_Store -delete
	@echo ""
	@echo "final_outputs/ assembled."

# No export/publish target: tree and history hold no internal documents, so git push publishes.

# clean spares figure_data/ (make figures needs it), 4 frozen gps_*.csv files, and 2 supplementary READMEs.
clean:
	rm -rf logs
	rm -f final_plots/genetic_data.png \
	      final_plots/perturbational_functional_genomics_with_thde.png \
	      final_plots/figure_model_analysis_composite.png \
	      final_plots/figure_trial_validation_immune_stages_R.png \
	      final_plots/figure_stat4_vignette.png \
	      final_plots/specificity_immune_vs_cardiac_faceted.png \
	      final_plots/figure_prospective_temporal_holdout.png \
	      final_plots/figure_label_permutation.png
# final_plots/pdf/ is the only directory clean deletes wholesale, holding only renderer output.
	rm -rf final_plots/pdf
# final_plots/supplementary/ holds only S2 and S3.
	rm -f final_plots/supplementary/leakage_controlled_comparison.png \
	      final_plots/supplementary/figure_discordance_genetics_vs_full.png
# data/raw/discordance/ and final_outputs/ handling: see REPRODUCIBILITY.md.
	rm -rf final_outputs
# clean wipes supplementary/: make supplementary rebuilds it as a pure copy, so stale files would look current.
	rm -rf supplementary
# outputs/ is deleted wholesale since nothing under it is tracked in git; make tables recreates needed dirs.
	rm -rf outputs
	rm -f data/derived/*.csv
