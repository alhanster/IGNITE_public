# Reproducibility notes

This document records what the build reproduces, what it does not, and every constraint that
affects whether a re-run yields the committed numbers. It is the companion to the Makefile:
the Makefile says how to run each stage, this says what the results mean and where they can
go wrong.

- [Environment and version pins](#environment-and-version-pins)
- [The two-stage build and the figure_data contract](#the-twostage-build-and-the-figuredata-contract)
- [Stage-1 step order](#stage1-step-order)
- [What the verification targets check](#what-the-verification-targets-check)
- [What a green gate does not cover](#what-a-green-gate-does-not-cover)
- [Opt-in targets](#optin-targets)
- [Statistical methods and seeds](#statistical-methods-and-seeds)
- [Data sources and frozen inputs](#data-sources-and-frozen-inputs)
- [Figure and table numbering](#figure-and-table-numbering)
- [Reproduction boundary](#reproduction-boundary)
- [Conventions](#conventions)

---

## Environment and version pins

### Python (stage 1, `make tables`)

Pinned versions in `requirements.txt` reproduce `make tables` and `make figures` byte-for-byte,
verified under Python 3.12.13 on macOS 26.5.1 arm64. Versions are pinned
exactly, not as ranges, because a minor version change can alter results:

- **xgboost 3.3.0** is the load-bearing pin. Under 3.3.0, `full_model_pu_scores.csv` reproduces
  to CSV precision across runs (max diff 1.1e-16, Spearman rho 1.0000000000). Under 3.2.0,
  results differ deterministically: STAT4 moves from rank 9 to rank 11 and panelD_recovery
  counts shift by plus or minus 1, because adjacent score gaps in the top 20 (0.00007 to 0.0027)
  are comparable to the 3.2.0-versus-3.3.0 score offset (about 2e-03). CV AUC (to three decimal
  places), the coalition ladder, and the Nadeau-Bengio stars are unaffected. xgboost 3.3.0 also
  publishes wheels only for Python 3.12 and later (3.0-3.2 supported 3.10+; 3.3.0 dropped 3.11
  support); on a Python 3.11 interpreter the pin cannot install, and pip reports the newest
  available wheel rather than the newest released one, so the failure can read as though 3.3.0
  was never released.
- **matplotlib** is pinned because it stamps its version into the PNG `Software` chunk: 3.11.1
  vs 3.11.0 differ in 5 bytes while drawing identical pixels.
- The `# python_requires: >=3.12` line in `requirements.txt` is parsed by both
  `tools/detect_python.sh` and `_version_guard.py`; the `>=X.Y` form must be preserved.
- An exact pin with no matching wheel for a given platform fails installation rather than
  silently installing a different version.

`src/analysis/common/_version_guard.py` runs `check_pins()` at import of `pu_target_model.py`, as step 1
of `make tables`; a mismatch aborts the run rather than emitting results that will not reproduce.
`IGNITE_ALLOW_VERSION_MISMATCH=1` downgrades this to a warning, for deliberate cross-version work
only; it must never be set for a run whose outputs are committed. `make check-versions` runs the
same guard standalone, alongside the R-side check.

`check_pins()` skips packages that are absent rather than failing on them, so meeting the Python
version floor is not sufficient to guarantee a working environment. `tools/detect_python.sh`
therefore selects an interpreter in two passes: it prefers one that already has the pinned
xgboost installed, falls back to one meeting only the version floor (so the guard reports the
package problem against a usable interpreter), and if none meets even the floor prints nothing,
letting the Makefile default to `python3` so the guard reports the version failure. It runs as a
standalone script rather than inline Makefile shell because `make` strips everything after a `#`
before parsing, and a `#` inside `$(shell ...)` produces an unterminated-call-to-function-shell
error.

### R (stage 2, `make figures`)

Required: **R 4.4.1**, parsed by `tools/check_r_versions.R` from a `#   R X.Y.Z` line in
`R-requirements.txt`. There is no import-time guard on the R side analogous to
`_version_guard.py`; a version mismatch produces a byte-different figure rather than an
install-time error. `tools/check_r_versions.R` runs from `make check-versions` before a build
(not at render time): exit code 0 means the environment matches the pins, exit code 1 means a
mismatch that could change rendered pixels, since figures are verified by byte-comparison and
ggplot2/ragg output differs across releases.

Exact package versions:

| Package | Version |
|---|---|
| ggrepel | 0.9.6 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |
| jsonlite | 2.0.0 |
| dplyr | 1.2.1 |
| ragg | 1.5.2 |
| scales | 1.4.0 |
| systemfonts | 1.2.1 |
| textshaping | 1.0.0 |
| gridtext | 0.1.6 |
| ggtext | 0.2.0 |

Install with exact versions:

```r
remotes::install_version("ggrepel", "0.9.6")
remotes::install_version("ggplot2", "4.0.3")
remotes::install_version("patchwork", "1.3.2")
remotes::install_version("jsonlite", "2.0.0")
remotes::install_version("dplyr", "1.2.1")
remotes::install_version("ragg", "1.5.2")
remotes::install_version("scales", "1.4.0")
remotes::install_version("systemfonts", "1.2.1")
remotes::install_version("textshaping", "1.0.0")
remotes::install_version("gridtext", "0.1.6")
remotes::install_version("ggtext", "0.2.0")
```

`ggtext`/`gridtext` back `fig_stat4_vignette.R` panel c, which uses `element_markdown()` /
`geom_richtext()` so a markdown-italic STAT4 can appear in a plain multi-line title (plotmath's
`atop()` applies an implicit script-style size reduction to nested arguments, producing uneven
title-line sizes). `gridtext` is pinned for the same reason as `systemfonts` and `textshaping`
are pinned for `ragg`: it directly determines glyph pixels rather than acting as an incidental
transitive dependency.

`tools/setup_env.sh` installs R packages in three passes. Passes 1 and 2 fetch only the current
CRAN release, so a pin older than the current release installs successfully but leaves the pin
unmet (for example, `systemfonts` and `textshaping` are pinned at 1.2.1 and 1.0.0 while CRAN
currently ships 1.3.2 and 1.0.5). Pass 3 fetches from the CRAN Archive via
`remotes::install_version` to close this gap; skipping it silently installs the wrong font stack,
which determines rendered output.

Two prerequisites are not installable from pip or CRAN and are checked before any install begins:
native headers (`ragg` requires freetype via `ft2build.h`, `textshaping` requires harfbuzz via
`hb-ft.h`) and a working C++ toolchain, since several R packages are C++ and one is
Objective-C++. CRAN's macOS binary packages bundle these libraries, so the gap only appears once
a pin forces a source build. The header probe uses the combined include set rather than one
directory per header, because `hb-ft.h` includes `ft2build.h` and probing harfbuzz alone can
report it missing on a machine where it is present; pkg-config is not assumed present, and
Homebrew include paths serve as the fallback. `IGNITE_CXX_SDK_FALLBACK=1` builds against the
SDK's libc++ when the default toolchain is broken; it is opt-in and offered only after a probe
confirms the SDK libc++ compiles.

On macOS, CRAN binary installs of `ggplot2`, `dplyr`, and `ragg` yield 4.0.2, 1.2.0, and 1.5.1
respectively, one patch behind the pins above; installing those three with `type = source`
avoids the mismatch. `ragg` additionally requires freetype, libpng, jpeg-turbo, and libtiff
headers, available via Homebrew. `gridtext` requires the CRAN `jpeg` package to build from
source, which requires libjpeg headers visible to the compiler; on macOS without pkg-config on
PATH, headers provided by jpeg-turbo (already required by `ragg`) may not be found, and
installing `jpeg` from a CRAN binary first resolves a dependency-not-available error.

System libraries used for glyph rasterization (freetype, libpng, libtiff, jpeg-turbo) cannot be
pinned by `R-requirements.txt`. The verified build linked freetype 2.14.3, libpng 1.6.58, libtiff
4.7.2, and jpeg-turbo 3.2.0 on macOS 26.5.1 (arm64). `ggrepel`'s label solver runs in compiled
C++ via Rcpp (1.0.14 in the verified build, itself unpinned), so panel a of the discordance
figure carries the same native-layer exposure as the PNG bytes. A different native stack can
change PNG bytes even when every pinned R package matches, so a byte mismatch does not
necessarily indicate a pipeline regression.

### Environment layout and overrides

Both stacks install into project-local, gitignored directories, `.venv/` and `.rlib/`, rather
than a shared or global environment, so the pinned versions of `ggplot2`, `dplyr`, and `xgboost`
do not collide with other projects on the same machine.

Stage 1 (`make tables`) requires Python; stage 2 (`make figures`) requires R. The two
environments may differ; override with:

```bash
make tables PY=/path/to/python3
make figures RSCRIPT=/path/to/Rscript
```

Both use `?=` in the Makefile, so an environment variable of the same name takes precedence over
the default. `PY` must point to the pinned environment (`xgboost==3.3.0`); with auto-detection in
use, `make check-versions` prints the resolved interpreter path rather than just its name, since
an auto-detected `python3` may resolve to an unintended conda base environment.

`make check-versions` validates both stacks before a rebuild: the Python check catches the
xgboost mismatch class also caught at import time, and is the only safeguard against a mismatched
ggplot2 or ragg version producing figure differences that resemble code regressions, since there
is no render-time guard on the R side.

---

## The two-stage build and the figure_data contract

`figure_data/` is the contract between the two build stages. `make tables` (stage 1) writes it;
`make figures` (stage 2) reads only from it and requires no model fitting, writing its rendered
output directly into `final_plots/` with no intervening publish step. Because the two stages
communicate solely through this directory, it is committed rather than gitignored, and a
`make tables` run leaves the working tree dirty under `figure_data/` — that is expected, not a
bug; `make verify-tables` is the signal to trust, not `git status`.

The contract is narrower than what any one script computes. `pu_target_model.py` exports, as its
Stage-2 file, only the two columns needed for the 60-bin histogram, in per-gene row order sorted
identically to `OUT_RANK`; the full ranking table stays unexported scratch under `outputs/`.
Similarly, `knn_p95` (the dashed reference line on panel d) is not a hand-maintained constant in
the seed meta but is read by `build_vignette_meta.py` directly from
`build_vignette_panelD_knn.py`'s output.

### Numbered supplementary tables

`tools/package_supplementary_tables.py` packages the numbered supplementary tables directly out
of `figure_data/`, never out of `final_plots/supplementary/`. This is sufficient for correctness
because `figure_data/` is exactly what `make verify-tables` checks, and it means `final_plots/`'s
sha256 baseline is untouched by table changes: adding or editing a supplementary table needs no
baseline refresh. These tables never enter `final_plots/` and so stay out of the sha256 baseline
by construction; reintroducing a publish step for them would also require updating
`verify-figures`, which currently treats a gitignored non-PNG output as a hard failure.

Table names and numbers live in `tools/supplementary_tables.tsv`, not in the packaging script,
because `final_plots/supplementary/` filenames alone cannot carry a numbering that shifts when a
table is added or removed. The packaging script validates that numbering for gaps and duplicates,
excludes the frozen GPS tables, and derives a column dictionary from each file's header. It runs
before `write_submission_readme.py`, which lists the directory it produces.

`supplementary_data/` itself is not assembled. The statistics it would have contained are
reported instead in Supplementary Tables 7, 9 and 11, sourced from `figure_data/`.

---

## Stage-1 step order

`make tables` runs 30 Python steps. **The order is load-bearing; the stage does not support
parallel execution with `-j`.**

Stage-1 scripts write to two places. `figure_data/` holds the committed contract checked by
`make verify-tables`. `outputs/` holds scratch: anything no downstream figure reads. Nothing under
`outputs/` is committed, and `make clean` removes it; a fresh clone has no cached gene or PU scores
and must run `make tables` to obtain them. `.gitignore` carries one negation per file that a
retained stage-1 script still writes into an otherwise-ignored directory; removing a negation
causes `make tables` to leave an untracked file in `git status` on every run.

The dependencies that make the order load-bearing:

- `reshape_attribution_to_panels.py` writes `panelA_significance.csv` and `panelC_marginal.csv`.
  `compute_nb_significance.py` must run after it and rewrites those same two files in place to
  install Nadeau-Bengio corrected p-values.
- `pu_target_model.py` writes `full_model_pu_scores.csv`, read by three later scripts:
  `build_vignette_meta.py`, `build_leakage_controlled.py`, and `build_ranked_atlas_table.py`.
  `build_leakage_controlled.py` also imports `delong_test` from `compute_delong_specificity.py`.
- `build_vignette_meta.py` reads `full_model_pu_scores.csv` (from `pu_target_model.py`, to keep
  STAT4's rank annotation current) and `vignette_knn_stats.json` (from
  `build_vignette_panelD_knn.py`, for the knn_p95 reference line in panel D). Both must run first;
  the script aborts rather than substituting a default value if either input is missing.
- `build_ranked_atlas_table.py` runs after `pu_target_model.py` and before its four downstream
  readers of `figure_data/ranked_atlas.csv`: `build_immune_stages_data.py`,
  `build_feature_dictionary.py`, `build_univariate_by_group.py`, and
  `build_orthogonality_tables.py`. Run out of order it fails on a clean tree, since
  `ranked_atlas.csv` does not yet exist. It performs a pure join of two already-committed inputs,
  `outputs/model/full_model_pu_scores.csv` and `data/perturbseq/pu/pu_model_matrix.parquet`, with
  no model fitting or randomness. Its output is not published into `final_plots/`; it is
  Supplementary Table 1, packaged from `figure_data/` by `tools/package_supplementary_tables.py`,
  which keeps it out of `final_plots.sha256`.
- `build_immune_stages_data.py` must follow both `build_ranked_atlas_table.py`, whose output it
  joins to add six reporting columns to Supplementary Table 10, and `make_panelD_recovery.py`,
  whose `genetics_only_pu_scores.csv` it reads for the genetics-only arm in panels b and c. The
  two arms share `pu_role` labels and differ only in features, so both rank the same 19,162
  non-approved genes with a different ordering.
- `build_feature_dictionary.py` must follow `pu_target_model.py`, for
  `pu_feature_importance_stability.csv`, and `build_ranked_atlas_table.py`, for
  `figure_data/ranked_atlas.csv`, which it uses to compute the observed per-feature correlation
  against the published ranking.
- `build_supp_stat_sheets.py` runs after `build_score_by_group.py`,
  `build_heldout_trial_auc.py`, `build_immune_stages_data.py`, `build_leakage_controlled.py`, and
  `build_vignette_panelD_knn.py`, whose JSON outputs it reshapes. It computes no new statistics:
  each value is copied from a committed JSON file, reshaped into a tidy frame, and checked against
  its source. Outputs are written to `figure_data/` rather than produced during packaging, so
  `make verify-tables` can check them against the committed values; the raw JSON files still ship
  via `final_plots/`.

  | Output | Rows | Table |
  |---|---|---|
  | `heldout_trial_stats.csv` | 9 | Supp Table 7, discrimination statistics |
  | `score_by_group_summary.csv` | 3 | Supp Table 7, per-group descriptives |
  | `score_by_group_tests.csv` | 6 | Supp Table 7, group comparisons |
  | `scrambled_control_stats.csv` | 106 | Supp Table 8, scrambled-feature control |
  | `immune_stages_bootstrap.csv` | 4 | Supp Table 10 / Fig 5b, both models x 2 trial-stage thresholds at depth 50 |
  | `immune_stages_depth_sweep.csv` | 6 | Supp Table 10 / Fig 5c, both models x 3 shortlist depths, Phase >= 1 only |
  | `leakage_controlled_delong_stats.csv` | 19 | Supp Table 11, paired DeLong per arm |
  | `vignette_knn_summary.csv` | 6 | Supp Table 14, neighborhood test |

  `heldout_trial_roc_curves.csv`, shown alongside `heldout_trial_stats.csv` in Table 7, is thinned
  to about 1,400 points per curve for plotting; a trapezoid integration over it yields 0.636481,
  matching the reported 0.636481 to machine precision. AUC and DeLong statistics come
  from `heldout_trial_stats.csv`, computed over the full curve, not from the thinned
  plotting sheet.

  `score_by_group_tests.csv` contains two Mann-Whitney rows for the same pair of groups: U =
  8,530,276.5 for the one-sided held-out test (in-trial first) and U = 4,871,968.5 for the
  two-sided pairwise test (non-target first). The two values sum to 727 x 18,435, and the
  p-values differ by exactly a factor of two. Holm correction applies only to the three pairwise
  tests; the headline held-out test, `kruskal_p`, and `anova_p` are raw, so `p_holm` is empty on
  those rows.

- `build_discordance.py` runs last. It asserts its refit reproduces the committed ranking
  (`spearman_vs_committed_full` exactly 1.0) and can abort on a stale STRING/GO cache.

`src/analysis/README.md` is the map for the PU model, which has the widest blast radius in this
stage.

---

## What the verification targets check

### `make check-versions`

`audit_coverage()` in `src/analysis/common/_version_guard.py` checks that every script constructing an
`XGBClassifier` calls `check_pins()`. An unguarded script could train under a different package
version and produce numbers inconsistent with the committed artifacts. It returns the list of
unguarded paths, empty when the repository is clean.

### `make verify-tables` and `make verify-figures`: the five-class system

Both targets sort every checked file into one of five classes before printing anything, then
report failures one class at a time (the baseline is sorted by path and classes are not
contiguous within it):

| Class | Definition | Severity | Why |
|---|---|---|---|
| MISSING | shasum cannot open the file | Fatal by default | Absence is not the same as drift. `VERIFY_ALLOW_MISSING` downgrades this to reported-only; only the `supplementary` target sets it. |
| RENDERED | a gitignored `.png` **or** `.pdf` | Soft, non-fatal | Non-fatal because a *different machine* can reproduce every number and still fail these hashes: `ragg` rasterizes through freetype, libpng and the installed fonts, which CRAN does not pin. On one machine under the pinned R stack PNG bytes **do** reproduce — 10 of 10 byte-identical across two independent clones — so a PNG mismatch here is a real signal, most often a render that did not resolve `.rlib` and used a different `ragg`. PDFs are in this class for a separate reason that no stack fixes: every R PDF device stamps a wall-clock `/CreationDate`. All 20 baseline files are in this class today. |
| DOCS | a `.md` file | Reported only, never fatal | Protected by git rather than by this baseline; covers the two supplementary READMEs. |
| PUBLISHED | gitignored, and not a `.png`, `.pdf` or `.md` | Fatal | A byte-deterministic build output copied into `final_plots/` from `figure_data/`, so any mismatch is real rather than a rendering artifact. No file currently falls into this class: the `figures` target renders only, and the numbered supplementary tables go straight from `figure_data/` to `final_outputs/` via `make final-outputs` without passing through `final_plots/`. It exists to protect any table added to `final_plots/` later. |
| COMMITTED | tracked, not a `.md` file, sourced externally with no producing script | Fatal | A mismatch means the file was corrupted, truncated, or overwritten. No file currently falls into this class; it exists to protect any file added later. |

PUBLISHED and COMMITTED are kept as separate classes because they fail for different reasons and
only one indicates corruption, and conflating them admits a byte-identical gitignored
build output to be reported as corruption.

Absence is tested before filename-based classification: shasum's open-failure and hash-mismatch
messages both contain `: FAILED`, and classification keys only on filename, so without an
explicit absence check first, a missing file would be silently classified as a
present-but-soft-mismatched RENDERED file.

`verify-tables` and `verify-figures` are kept as separate pass/fail results rather than combined.
`verify-tables` confirms the numbers reproduce (hard failure, portable across machines).
`verify-figures` confirms the pixels reproduce (reported only, not gating), because gating on
pixel bytes would fail the build on any rendering difference while leaving `figure_data/` - which
holds every AUC, p-value, star, rank and count reported in the manuscript - unverified. Keeping
the checks separate lets a reviewer on a different system reproduce every reported number exactly
even when pixel output differs.

### PDF outputs and the supplementary target

`final_plots/pdf/` holds only generated vector output for Figs 2-8 and S1-S3 and is covered by a
single gitignore directory rule rather than per-file entries, since the directory can hold no
tracked files. It is the journal's vector deliverable, written by the same renderers as the
adjacent PNGs, and is ignored for that reason plus one more: each PDF embeds a wall-clock creation
timestamp, so it differs on every render. `verify-figures` classifies these files as RENDERED
(soft, non-exact) for that reason.

`vector_output.R` writes PDFs only through the quartz device, available on macOS. On any other
host each renderer writes its PNG, logs that the PDF was skipped, and carries on, so `make
submission` completes on Linux (including WSL) with every table and number intact. Off macOS,
`check-versions` reports a missing Arial Unicode MS font rather than failing: PNG pixels then
differ from the macOS baseline, which `verify-figures` already treats as a soft difference. Setting
`IGNITE_SKIP_PDF=1` forces the same path on macOS, which is how it is tested. A skipped render
also deletes any older PDF of that figure, so a PNG never ships beside a stale PDF.

Absent PDFs are all-or-nothing. When `final_plots/pdf/` holds no PDF at all, `verify-figures`
reports the missing baseline PDFs as not rendered rather than failing, `final-outputs` packages
no `pdf/` directories, and the generated `final_outputs/README.md` says so. A partial set of PDFs
indicates an interrupted or failed render: `verify-figures` treats the gaps as missing baseline
files, and the supplementary packager aborts.

The `supplementary` target renders and packages only the supplementary display items rather than
running the full `figures` and `final-outputs` builds. It writes to its own `supplementary/` tree
rather than `final_outputs/`, since `final-outputs` removes and rebuilds `final_outputs/` from
scratch and sharing a tree would make build ordering matter. Verification runs inside this recipe
after rendering, not as a prerequisite, since a prerequisite would compare against stale pixels
from before rendering. The recipe sets `VERIFY_ALLOW_MISSING=1` (no other target sets this
variable): it renders three of the ten figures, so the other seven are reported as missing rather
than failing the build. The relaxation changes only the exit code; the listing of missing figures
is unaffected.

### Cross-platform reproduction

All results reported in the manuscript and supplementary tables were generated on macOS (Apple Silicon) using the environment provided in this repository. Minor numerical differences may occur when the analyses are reproduced on other operating systems or hardware, because of differences in numerical libraries, floating-point operations and parallel computation. These differences do not affect the qualitative conclusions or statistical interpretation of the analyses. In the exploratory discordance analysis (Supplementary Fig. S2 and Supplementary Table 9), the core gene set and features near the FDR threshold may differ.

The committed `figure_data/` tables, and every number in the manuscript, were produced on Apple
Silicon macOS (arm64). That is the reference platform: there, `make tables` reproduces all 98
tables byte-identically, and `verify-tables` fails hard on any difference.

On any other platform the model-derived tables drift. XGBoost and the linear-algebra libraries it
calls produce slightly different floating-point results on Linux x86 than on Apple Silicon, and
PU bagging fits thousands of small models, which amplifies those differences. Version pins cannot
remove this. Tables that do not depend on the model differ, if at all, only in the last digit of a
printed float. Off the reference platform, therefore:

- `verify-tables` reports which tables differ and exits successfully instead of failing.
- Checks that tie a table to the committed run print a warning instead of stopping the build:
  the controlled-arm values pinned in `build_leakage_controlled.py`; the permutation-null tables
  from `make permutation-null`, which `build_supp_stat_sheets.py` checks against the refitted AUC
  ladder; and the STRING/GO tables from `make discordance-network`, which `build_discordance.py`
  checks against the current core gene sets. Those opt-in tables are not rebuilt by `make
  tables`, so off the reference platform they still describe the reference run.
  `src/analysis/common/_platform.py` holds the one platform test these checks share.
- The generated `final_outputs/README.md` states the platform it was built on.

#### Measured drift on Linux

Measured on 2026-09-12 with a full `make setup` and `make submission` on GitHub Actions: Ubuntu 24.04
x86_64, Python 3.12.14, R 4.4.1. Both commands completed, all 40 steps.

54 of the 98 tables differ from the committed reference. In 10 of them the only change is the last
digit of a printed float (under 1e-9). The other 44 carry model-derived values.

**No significance star or label changed in any table, and no starred result in any figure crossed
a threshold.** Among the unstarred values:

- **S2 panel c draws 10 features rather than 9.** The panel shows features at FDR < 0.05, and on
  Linux the demoted `n_sig_regulated_cytokines_Stim8hr` moves from 0.057 to 0.047. This follows
  from the S2 core set, which is 49 genes on the reference platform and 52 on Linux (AOAH leaves;
  EGFL6, HMGXB4, HYI and SS18L1 join). The core set is defined by quadrant cuts at one seed, and 10
  of its 49 reference genes sit within 0.01 of a cut, so a change of this size is expected.
- **Other raw p-values that cross 0.05, 0.01 or 1e-3 appear only in supplementary sheets:** the
  remaining S2 feature and stratified-validation tests (Supp Table 9); the genetics-only Phase ≥ 3
  enrichment at depth 100 (Supp Table 10, not plotted; its count goes from 6 to 5); and the
  one-sided T0 = 2013 DeLong test (Supp Table 12).
- **All 10 PNGs differ from the reference baseline**, as expected with different fonts and
  rasterization libraries. Linux produces no PDFs.
- **Tables from the opt-in targets are not rebuilt** by `make tables`, so a Linux package ships
  them as committed: the permutation null and scrambled-feature control (Supp Table 8 and S1) and
  the STRING/GO network tables (S2).

| Quantity | Reference (Apple Silicon macOS) | Linux x86_64 |
|---|---|---|
| CV ROC-AUC, full model | 0.7876 | 0.7877 |
| CV ROC-AUC, genetics only | 0.7682 | 0.7681 |
| Fig 4b Shapley share, genetic / perturbational / observational (%) | 68.1 / 12.6 / 19.3 | 68.1 / 12.5 / 19.4 |
| Held-out in-trial AUC, full model | 0.6365 | 0.6367 |
| Held-out in-trial paired DeLong p | 0.0122 | 0.0116 |
| Fig 5b BH p, full model, Phase ≥ I | 0.0004 | 0.0004 |
| Fig 5b BH p, genetics only, Phase ≥ I | 0.0397 | 0.0386 |
| T0 = 2014 ΔAUC / two-sided DeLong p | 0.0364 / 0.0002 | 0.0361 / 0.0002 |
| T0 = 2015 ΔAUC / two-sided DeLong p | 0.0495 / 0.0454 | 0.0488 / 0.0472 |
| Fig 6b emergent genes the full model ranks higher, of 134 | 94 | 93 |
| STAT4 score / rank | 0.9633 / 13 | 0.9643 / 11 |
| S3 leakage-controlled AUC / DeLong p | 0.6451 / 0.0003 | 0.6453 / 0.0003 |
| S2 core genes, promoted / demoted | 36 / 13 | 37 / 15 |
| S2 Spearman rho between arms | 0.832 | 0.833 |
| S2 panel c features at FDR < 0.05 | 9 | 10 |

### Build-time assertions inside analysis scripts

Several scripts run under `make tables` carry their own pass/fail checks, independent of the
`verify-tables`/`verify-figures` baseline comparison:

- **`build_discordance.py`**: the full arm (genetics-only plus a 19-feature functional-genomics
  block, fit at `pum.SEED` and `T=pum.T_BAG`) is asserted to reproduce
  `outputs/model/full_model_pu_scores.csv` exactly (`spearman_vs_committed_full == 1.0`); the run
  aborts otherwise, since seed, `T_BAG`, tree depth and the xgboost pin are all certified by this
  one check. FDR correction uses `scipy.stats.false_discovery_control(method="bh")` in place of
  `statsmodels.multipletests(method="fdr_bh")`; the two agree to within 6e-16 on the source
  p-values with no gene crossing the FDR 0.05 or 0.10 threshold differently. HGNC family
  enrichment is not computed, since it is not a valid statistic for this gene set. The
  GO-dark/UniProt dark-proteome comparison is computed in `fetch_discordance_network.py` instead,
  avoiding a dependency on a gitignored file in this script's verification path.
  The staleness gate runs before the network fetch step: if a cache is present but describes a
  different gene set than the current core sets, the run aborts (this is the only point where
  drift between the cache and the current model configuration is detectable); if the cache is
  absent, the `net_*`, `go_*` and `dark_*` keys are simply omitted, and the renderer's own gate
  stops the corresponding panels. Each of the four cache checks (network, GO, dark-proteome,
  matching input files) is gated independently at module level by its own input file's existence,
  so a partial checkout cannot skip validation while still folding stale scalars into `discordance_stats.json`.
  The GO cache has no independent staleness check; it is validated by the invariant that the
  number of foreground genes with a GO biological-process term plus the number counted as GO-dark
  equals the promoted set size (28 + 8 = 36), computed off the same GAF as both the GO and
  dark-proteome tables. A missing dark-proteome table raises an error rather than skipping the
  check, since without it nothing pins the GO cache to a gene set.
- **`build_drug_status_table.py`**: the per-gene rollup is an independent flattening of the same
  evidence used to build the per-drug table. Agreement with the rollup on its two untruncated
  columns, `n_drugs` and `furthest_stage`, checks that the multi-target explode step and the stage
  ordering are correct.
- **`build_heldout_trial_auc.py`**: class counts for held-out in-trial genes (727) and non-target
  genes (18,435) are computed at build time and written, not hardcoded in prose. The genetics-only
  comparator is re-scored at `pum.SEED` and `pum.T_BAG` using config shared by import with the
  model-analysis recovery panel, and the build asserts this genetics-only ranking reproduces the
  recovery counts in `panelD_recovery.csv`. The check imports its counting helpers from
  `make_panelD_recovery` rather than reimplementing them, so it verifies agreement across two
  independent re-derivations rather than a duplicated counting routine agreeing with itself.
- **`build_immune_stages_data.py`**: this script and the atlas both derive scores from
  `full_model_pu_scores.csv`, so the two must be bit-identical; a mismatch indicates one side read
  a stale copy of the file.
- **`build_vignette_meta.py`**: `knn_p95` and the panel c GWAS counts are required inputs and fail
  loudly if missing rather than falling back to a substitute value, since their absence indicates
  a broken checkout rather than an unrun opt-in fetch, and `verify-tables` byte-compares this file
  against the committed output.
- **`build_supp_stat_sheets.py`**: Table 10 packages the top-50 nominations sheet and its
  bootstrap sheet in one workbook; the caption's claim that nine of the 50 carry an
  immune-indicated drug in trial is checked by confirming the two sheets agree on that count and
  that the observed folds re-derive. The comparison sheet rounds values to three decimal places
  while the DeLong JSON does not, so its check applies a rounding tolerance rather than exact
  equality - this is what would catch the two sheets describing different runs, the failure mode
  `PROVENANCE.md` records for the frozen pre-depth=2 GPS tables. The renderer draws the
  significance star string verbatim, so the star is re-derived from `p_holm` rather than raw `p`
  (the two arms form one correction family and the marker is decided on the corrected value); both
  give the same star here, so this check verifies only that the referenced column matches the one
  the producer used, not that the threshold itself is correct. For Table 14, the condition
  assertion is load-bearing: `data/perturbseq/knn/` must hold exactly one signature file,
  `signatures_Stim48hr.npy`, matching the caption's claim that the signatures come from that
  condition; a future pooled-condition rerun must update this check explicitly rather than
  silently changing whether the caption's claim is true.

---

## What a green gate does not cover

`make tables` rewrites 85 of the 96 tables in `figure_data/`. The remaining 11 belong to opt-in
targets that `make tables` does not run: the 6 `discordance_*` caches, the 3 permutation-null
files, and the 2 `vignette_panelC_gwas_*` files. `verify-tables` passes these trivially on an
ordinary clone, since nothing rewrote them to compare against the commit; a pass proves only that
they are unchanged, not that they are correct.

**Permutation null.** `make tables` does not rewrite the `label_permutation_*` files, so
`verify-tables` again shows only that they have not changed. They were, however, regenerated
under the current pins on 2026-09-06 against the 24-feature matrix — 1,000 permutations plus 99
scrambled-feature draws, 8 h 33 m wall-clock at `PERM_JOBS=11`. The committed draws are
therefore a product of this toolchain, not an inheritance from an older one. Table 8's two
pass-through sheets are checked against each other and against `panelA_auc_ladder.csv` by
`build_supp_stat_sheets.py`: all four `emp_p` values re-derive from the committed 1,000 draws as
`(k+1)/(n+1)`, and the four `observed` levels match. What remains uncovered is repetition: the
null has been generated once under these pins, not shown to re-derive. Only a second
`make permutation-null` would close that.

**KNN perturb-seq signatures.** The signature matrix
(`data/perturbseq/knn/signatures_Stim48hr.npy`, 294 MB) is gitignored and absent on a fresh clone.
`build_vignette_panelD_knn.py` treats that absence as a skip rather than an error, so `make tables`
and `make verify-tables` complete offline, and `vignette_panelD_knn.csv` /
`vignette_knn_stats.json` stay at their committed values. On such a clone `verify-tables` passes
both files trivially. Only `make knn-signatures` followed by a rerun of the script re-derives
them; the committed row and column gene-name lists let their identity be checked without the
matrices themselves.

**Discordance network cache.** `validate_against_cache` recomputes the induced subgraph for the
gene sets in `discordance_network_nodes.csv` and compares edge count, connected-node count, and
the exact edge list against the committed `discordance_network.csv` and
`discordance_network_edges.csv`. This can only confirm agreement with the cache it is itself
built from; it cannot detect a STRING release with different scores or identifier mapping than
the one used to build the cache, which is exactly when the cached p-values stop reflecting current
STRING data. Separately, the edge cache filename encodes only `MIN_SCORE`, not the gene universe:
if the gene universe changes without a matching cache rebuild, `load_edges()` silently returns a
graph missing edges that touch the new genes, and `validate_against_cache` cannot see this either,
since it recomputes from that same cache. `make clean` preserves `data/raw/discordance/`, so this
stale state survives a clean.

**Final package.** `make verify` checks `figure_data/` and `final_plots/`, not `final_outputs/`.
`tools/compare_submission.py` rebuilds the `final_outputs/` package from a second clone and
compares it against the package being submitted, covering that gap. Two of the four output
formats embed a wall-clock timestamp at write time and so differ on every rebuild regardless of
content: `.xlsx` (`docProps/core.xml`, `dcterms:created`/`dcterms:modified`) and `.pdf`
(`/CreationDate`, `/ModDate`, and the `/ID` derived from them). Everything else in those two
formats is compared by SHA-256 after the timestamp fields are blanked. `final_outputs/README.md`
lists a size and SHA-256 for every file; its rows for the `.xlsx` and `.pdf` files are normalised
before comparison, but its prose and its `.png`/`.csv` rows are compared exactly, so a content
change there still surfaces. `.png`, `.csv`, and `.json` outputs carry no exemption. The script
exits 0 only if every file in the package reproduces under these rules.

**Baseline provenance.** `tools/write_baseline.sh` reuses the existing version checkers rather
than duplicating them, and both exit non-zero on mismatch. But the versions stamped into
`final_plots.sha256` describe the machine that ran the script, not necessarily the one that
rendered the figures being hashed. A baseline regenerated on a mismatched machine over an
otherwise correct tree records correct hashes with incorrect provenance, and is silently
authoritative from then on unless its header states the mismatch.

**Files outside the manifest.** `shasum -c` checks only files listed in the baseline manifest, so
any file added to `final_plots/` afterward, including a fabricated table, is invisible to
verification and can reach `final_outputs/` unchecked. An unlisted `.md` file is treated as a soft
or expected draft, but only because git preserves a copy of it; an unlisted file that is neither
hashed nor tracked has no recovery path and must not be reported as protected. Unlisted data or
figure files are reported as unvouched deliverables.

**Group counts in `build_score_by_group.py`.** The output CSV is filtered to genes with a
cross-fitted `rank_pctile_cv`, but the statistical tests run on the unfiltered frame, and the group
counts written to the JSON (panel d's axis labels) come from that unfiltered frame. The script now
asserts that `rank_pctile` and `rank_pctile_cv` agree on which genes are scored, so a divergence
that would overstate the plotted n fails the build rather than passing quietly; upstream,
`pu_target_model.py` also asserts that the cross-fit scores all 19,502 rows.

---

## Opt-in targets

Six targets sit outside `make tables`, which must stay offline and complete in minutes. Each is
invoked explicitly:

```bash
make permutation-null            # ~2-3 h, ~3M XGBoost fits; checkpoints and resumes
make permutation-null-finalize   # seconds; recompute p-values from the existing null
make discordance-network         # ~10 min first run + a one-time ~135 MB download
make vignette-panelc             # seconds, network; re-derives STAT4 panel c from Open Targets
make gps-inputs                  # seconds, network; re-derives the GPS inputs from the published table
make knn-signatures               # ~3-4 min, network; streams the perturb-seq knockdown signature matrices
```

`permutation-null`, `discordance-network`, `vignette-panelc` and `gps-inputs` write straight into
the committed `figure_data/` contract and have no dry-run mode. A smoke run
(`make permutation-null PERM_N=8 PERM_SCR=2`, `make discordance-network DISC_NULL=5`) leaves
partial output in the tree; revert with `git checkout -- figure_data`. As with `make tables`, a
a real run leaves `figure_data/` dirty; `make verify-tables` is the check that the rewrite
reproduced, not `git status`.

Three optional Python dependencies are unpinned by policy because no committed result depends on
their exact version, and each is imported only by one of these targets: `h5py`
(`build_knn_signatures.py`), `requests` (`fetch_vignette_panelC_gwas.py`), `openpyxl`
(`build_gps_per_gene.py`). A failed install of any of the three does not block reproducing
`make tables`. `pytest` is unpinned for the same reason: it backs `make test` only, is imported
by no build step, and its absence blocks neither `make tables` nor `make figures`.

### permutation-null / permutation-null-finalize

1,000 permutations by 4 rungs by 5 seeds by 5 folds by 30 bags is about 3 million XGBoost fits,
measured at 18.5 seconds per permutation at `--jobs 4`; budget 2-3 hours at `PERM_JOBS=11`. The
run checkpoints every 25 permutations into gitignored `outputs/permutation/`, so an interrupted
run resumes rather than restarting. The observed ladder is read from the committed
`figure_data/panelA_auc_ladder.csv`, not the gitignored `outputs/` copy, so the target runs from a
fresh clone and the null is calibrated against the same ladder draws as panel A.
`permutation-null-finalize` only recomputes p-values from the existing checkpointed null and takes
seconds.

### discordance-network

Builds Supplementary Table S2's STRING coherence, GO over-representation and dark-proteome
statistics. It reads about 135 MB of external archives: STRING v12.0 protein links (83 MB) plus
STRING's protein.info and protein.aliases maps, the UniProt-GOA GAF (11 MB), and a UniProt
annotation-depth table (3.5 MB). These archives are not committed; the eight small tables they
produce are, and they are the only `figure_data/` entries `make tables` does not rewrite.
Downloads land in gitignored `data/raw/discordance/` and are reused on later runs (`make clean`
does not remove them): about 10 minutes on first run, about 1 minute after. The target needs no
prior `make tables` run, since it reads the already-committed `discordance_scatter.csv` for the
`cov_class` and `gen_dec` columns used in covariate matching, and so runs from a fresh clone.

`make verify-tables` confirms the eight files are byte-identical to what is committed but cannot
confirm they still re-derive from the raw archives; this limitation is recorded in
`figure_data/PROVENANCE.md`.

The STRING null model matches `--n-null` gene sets to the core set on functional-genomics
coverage class by genetics-only decile, drawn from non-core genes, because STRING connectivity and
inclusion in the perturb-seq data both correlate with how well-studied a gene is. The reported
p-value is the one-sided empirical fraction of null sets with at least as many edges, with a +1
correction; at 500 draws the floor is p = 0.002. GO over-representation uses as background every
non-discordant gene measured in both functional-genomics blocks, so there is no background
sampling error; two backgrounds are written to `discordance_go_enrichment_full.csv`:
`full_universe` (all annotated non-discordant genes, approximately 16,961) and
`coverage_restricted` (approximately 8,059 annotated both-coverage genes, the primary background).

Rerunning is needed only when the core gene sets move, for example from a change to `pum.SEED`,
tree depth, or the feature matrix. `build_discordance.py` detects this and aborts with the
recipe: `make tables` (aborts, having written the new core sets), `make discordance-network`,
`make tables` again (passes, folding the new tables into `discordance_stats.json`). The committed network/GO
cache file is read, not regenerated, by `build_discordance.py`. `--validate` is the pre-flight
check: it recomputes the induced subgraph for the gene sets in the committed cache and confirms it
still reproduces the committed edge counts, which is what catches an upstream STRING release whose
scores or identifier mapping have moved, the case where overwriting the cache would silently
change the published p-value.

**What S2 does not show.** The quadrant sets are not target-enriched: one held-out trial gene
appears across both directions. The validation claim concerns shift direction within genetic
strata, which is a different statement, and quoting the marked genes as a target-enriched set
overstates it. This caveat lived in the repository README until 2026-09-12 and was moved here
when that section was retired.

**The GO background is a full set, not a sample.** The primary background for Supplementary
Table 9's GO sheet (the figure no longer draws GO) is the 8,059 coverage-restricted annotated
genes, the whole set rather than a draw from it, so no background sampling error applies. The
foreground `n_fg` is 28, not 36: eight of the 36 promoted genes carry no BP term, which the
dark-proteome test measures directly (22.2% versus 9.5%, OR 2.71, p = 0.018). The unrestricted
16,959-gene variant is not a uniform inflation of the restricted one, so the two backgrounds are
not interchangeable within a single statement.

**`has_cytokine` is a measurement flag.** In S2 panel c it records whether cytokine output was
measured for a gene, not a biological property of the gene, and should not be read as one.

These two notes also moved here from the repository README on 2026-09-12.

### vignette-panelc

Re-derives the STAT4 panel c data from Open Targets via `fetch_vignette_panelC_gwas.py`, the only
live-API query among these targets. It asserts Open Targets release 26.06 and aborts under any
other release, since association scores move between releases and this file backs a published
panel; that abort is the correct outcome, not something to work around.

### gps-inputs

Regenerates `data/specificity/gps_per_gene.csv` and `gps_drug_labeled_genes.csv` from the
published GPS supplementary table via `build_gps_per_gene.py`; these files have no other producer
in the repository. The source, `data/raw/gps/GPS_allgenes.xlsx` (about 5.6 MB), is gitignored
under `/data/raw/*`, so a clean clone lacks it and the target requires network access; this is why
it is excluded from `make tables`, which must build offline. The two outputs are committed, so
`make tables` reads them directly without touching the xlsx. The outputs feed
`build_leakage_controlled.py`. After running `gps-inputs`, `make tables && make verify-tables`
must reproduce `gps_max_overall` and the 453-gene list byte-identically; a mismatch indicates the
source GPS table or a derivation rule has changed, and every number in section 4.5a needs
re-checking.

### knn-signatures

Streams the perturb-seq knockdown signature matrices from the primary source,
`s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`, the genome-scale
CRISPRi Perturb-seq resource in primary human CD4+ T cells of Zhu, Dann et al. (2025), released
via the CZI Virtual Cells Platform. The file holds 33,983 perturbation-by-condition rows across
10,282 measured genes. Only the zscore layer (log2FC / lfcSE) is read; the layer must be
contiguous and uncompressed, checked by an assertion before reading. For each culture condition,
de-duplicated rows are kept where `ontarget_significant` is true and `n_cells_target >= 50`.

Outputs, written to `data/perturbseq/knn/`, are `signatures_<COND>.npy` (float32,
n_perturbed x 10282), `signatures_<COND>.genes.txt` (row gene symbols), and
`var_gene_names.txt` (column gene symbols). All three conditions are written and all three are
consumed downstream: `build_vignette_panelD_knn.py` reads `Stim48hr` alone for panel d of the
STAT4 vignette, and `build_knn_nearest_target.py` reads all three for Supplementary Table 15. The
`.genes.txt` files and `var_gene_names.txt` (110 KB total) are committed and overwritten by this
target; a clean `git status -- data/perturbseq/knn` after a run confirms the source has not
moved. The matrices themselves are not committed: `signatures_Stim48hr.npy` alone is 294 MB, and
serving it through Git LFS would consume a large share of a 1 GB/month free bandwidth budget, so
it is regenerated on demand instead.

`build_vignette_panelD_knn.py` has committed outputs (`vignette_panelD_knn.csv`,
`vignette_knn_stats.json`) and skips with a message when the signature matrices are absent, so a
clean clone reproduces `make tables` without ever running `build_knn_signatures.py`. On such a
clone, `make verify-tables` passes for those two files, but the pass is trivial: nothing rewrote
them, so it does not confirm real reproduction. Running `knn-signatures` followed by the panel-d
build step re-derives them for real.

### supplementary

`make supplementary` renders and packages only the supplementary display items, writing to its
own `supplementary/` tree rather than into `final_outputs/`, since `make final-outputs` removes
and rebuilds that tree and the two targets' write order would otherwise matter. Verification runs
inside the recipe, after rendering, rather than as a prerequisite, so that it checks the newly
rendered pixels rather than a previous build's output. The target renders three of the ten
figures and sets `VERIFY_ALLOW_MISSING=1`, so the seven figures it does not render are reported as
missing rather than failing the build.

---

## Statistical methods and seeds

### Cross-validation and seeds

The Shapley attribution decomposition (`attribution_decomposition.py`) evaluates all 8 coalitions
of the GENETIC, PERTURBATIONAL, and OBSERVATIONAL feature groups under 5-fold PU-bagging cross-validation,
averaged over 5 seeds. The empty coalition is chance (AUC 0.5); the exact Shapley value of each
group is computed from the 8 coalition scores, and the shares sum to the full model's
AUC-above-chance. Panels a, b, and c of Figure 4 report the 5-seed mean for aggregate metrics that
name no individual gene. `SEED=4` alone drives panel e recovery and every score-driven figure
(STAT4, trial-validation, specificity), because these require one fixed gene ranking rather than
an average.

`recompute_block_signal.py` (panel f) also reports a 5-seed mean, matching panels a-c; `SEED` there
fixes only the prior-score median split. It runs with `T_BAG = 60` rather than production's 200,
since lifts are stable well below 200 bags across the 7-model by 5-seed configuration.

In `coalition_auc.csv`, the chance row's `sd` of 0.0 is a definitional placeholder, not a computed
dispersion: `cv_auc` short-circuits to five copies of 0.5 per seed for the empty feature list, so
the 25-sample vector is 25 copies of 0.5. Every other `sd` is a population standard deviation
(ddof=0) over 5 folds by 5 seeds.

### Nadeau-Bengio corrected significance (Figure 4 panels a and c)

Panels a and c compare nested feature coalitions by cross-validated AUC, paired across folds.
Pooling the r*k fold scores into a naive paired test understates variance because folds share
training data (Dietterich 1998; Nadeau and Bengio 2003). `compute_nb_significance.py` replaces the
naive Wilcoxon signed-rank, fold-pooled p-values written by `reshape_attribution_to_panels.py`
with a Nadeau-Bengio corrected resampled t-test over r=5 by k=5 resampling:

```
corrected_var = (1/n + rho/(1-rho)) * sample_var(fold_diffs), with n = r*k = 25 and rho = 1/k
t = mean(diff) / sqrt(corrected_var), df = n - 1
```

The script reuses `attribution_decomposition`'s `cv_auc`, feature groups, and base learner, so the
refit scores identically to the bars it corrects. Point estimates are unchanged (deltas remain the
5-seed means from the reshape step); only significance values and panel c's standard error change.
It runs after `reshape_attribution_to_panels.py`, reading and rewriting `panelA_significance.csv`
and `panelC_marginal.csv` in place. The full-versus-genetic comparison moves from 1.1e-06
(naive Wilcoxon) to 8.1e-03 (Nadeau-Bengio); the uncorrected Wilcoxon values must not be published.
No published result reads a naive fold-pool p directly.

### Multiple-testing correction

Correction is applied within each pre-specified family rather than across the study.
`holm()` corrects small sets of pre-specified comparisons where each must individually hold;
`bh()` corrects screens, where the claim concerns the set as a whole, wrapping scipy's
`false_discovery_control` so every FDR calculation in the repository shares one implementation.
Both functions return a list positionally aligned with the input, since callers write the adjusted
value back into the same row as the raw p-value.

| Call site | Family | Method | Column |
|---|---|---|---|
| `make_enrichment_table.py` | 4 group x background IEI Fisher cells | Holm | `p_holm` |
| `compute_delong_specificity.py` | 3 pairwise ROC tests, per panel | Holm | `p_holm` |
| `build_leakage_controlled.py` | controlled and standard arms (2) | Holm | `p_holm` |
| `build_score_by_group.py` | 3 pairwise drug-status Mann-Whitney tests | Holm | `p_holm` |
| `build_univariate_by_group.py` | 3 Figure 2b Fisher tests of carrying any immune-GWAS signal | Holm | `p_holm` |
| `build_univariate_by_group.py` | 3 pairwise Figure 2b Mann-Whitney tests, genes with signal only | Holm | `p_holm` |
| `build_univariate_by_group.py` | 3 pairwise Figure 2c missense-Z Mann-Whitney tests | Holm | `p_holm` |
| `build_univariate_by_group.py` | 4 Figure 3a subset Fisher cells, per target class (2 families) | Holm | `p_holm` |
| `build_univariate_by_group.py` | 3 pairwise Figure 3b drug-status Mann-Whitney tests | Holm | `p_holm` |
| `build_univariate_by_group.py` | 2 Figure 3c/3d target-class Fisher tests, per condition (6 families) | Holm | `p_holm` |
| `build_immune_stages_data.py` | 8 Figure 5 enrichment tests | BH | `*_bh` |

`build_univariate_by_group.py` defines six kinds of family, twelve families in all. The first,
the Figure 3a row, is the one place here whose eight tests are corrected as **two** families
rather than one. The approved and in-trial classes are pre-specified as separate families of
four subsets each, because the Figure 3a claim is made per class — approved targets
are enriched in Th1, in-trial targets are enriched in Th1 — and not as a screen over the
(subset x class) grid, which is the condition the table above separates `holm()` from `bh()` on.
Treating the eight as one family instead is more conservative and changes no conclusion in this
panel: Th1 remains the only significant subset for both classes either way, and the two nearest
cells move only from `p_holm` 0.080 to 0.160. The alternative is recorded here because the
choice is visible in a published number and a reader recomputing per-subset, with no correction
at all, would reach a different conclusion — two cells sit at raw p ≈ 0.027.

Its second family definition, the Figure 3b pairwise row, is deliberately a near-duplicate of
`build_score_by_group.py`'s: both are three pairwise drug-status Mann-Whitney tests corrected
together, one on the residualized trans-regulatory burden and one on the model score. They are
separate families because they are separate claims on separate quantities, and both rows appear
above so that neither looks like the other's typo.

Its third, the Figure 3c/3d Fisher row, stops at the condition boundary for the same reason its
Figure 3a row stops at the target class: the panels facet by condition, and each condition is a
separately pre-specified question rather than one rung of a single claim. That boundary is the
one correction choice in this repo that **changes a conclusion** — the sole significant cell,
3d/Stim8hr/approved, is `p_holm` = 0.030 within condition and 0.183 under a twelve-test family,
where nothing survives. It is recorded in `figure_data/PROVENANCE.md` alongside the sheet, with
the note that the boundary followed the facet structure and predates the tests.

The remaining three kinds cover Figure 2. The Figure 2c row is an unremarkable sibling of the
Figure 3b one — three pairwise Mann-Whitney tests on a continuous feature, corrected together.
Figure 2b is the interesting case: its two rows test the **same feature** and are still separate
families. `gwas_score` is fully observed and never imputed, so a zero is a measured absence of
association rather than a gap; the panel drops those genes from its violin, and the fraction
dropped is itself a result. So the Fisher row tests presence on all genes and the Mann-Whitney
row tests magnitude on the genes with signal — different quantities, hence different families.
Correcting all six together was checked and changes nothing: five of six stay significant and no
cell crosses 0.05. The two blocks are worth reading side by side rather than separately, because
they disagree in a way that matters: approved and in-trial are indistinguishable in *whether*
they carry signal, but differ in *how much*.

Twelve families in one module is more than any other call site here, and the concentration is
structural rather than accidental: this module reads only the committed per-panel sheets of
Figures 2 and 3, so a test placed here is computed from the values the panel plots and cannot
drift from the figure. Every panel of both figures now has its statistics here except Figure 2a,
whose Fisher enrichment `make_enrichment_table.py` writes upstream because it needs the IEI gene
list rather than a plotted column. The families are not folded together because they test six
different quantities across two figures, and a single Holm across all of them would correct a
Figure 2b presence test against a Figure 3b rank test.

Left uncorrected:

- `compute_nb_significance.py`: the three coalition rungs are nested comparisons against a shared
  baseline, not independent hypotheses, and the Nadeau-Bengio correction already accounts for fold
  dependence.
- `build_temporal_holdout.py`: the T0 freeze-year sweep is a robustness check on a claim reported
  elsewhere.
- `_label_permutation.py`: p-values come from permutation nulls; three of the four level p-values
  sit at the 1/1001 resolution floor.
- `build_orthogonality_tables.py`: the 95-test grid supports a claim of absence of association,
  which correction could only strengthen.
- `build_discordance.py`: the four per-stratum Fisher tests are pooled by Mantel-Haenszel instead;
  the pooled p-value (`mh_p`) is the inferential statistic, and the per-stratum p-values are
  descriptive.

`build_discordance.py` and `fetch_discordance_network.py` call scipy's `false_discovery_control`
directly rather than through `bh()`; their `fdr` columns are checked against fixed values by
`verify-tables`.

For the enrichment table's star ladder, stars apply to `p_holm`; the four group x background
cells form one family. Holm leaves the four stars unchanged (it multiplies the largest p by 1),
which is why Figure 2 does not move under it. A Bonferroni correction would move it: multiplying
by 4 takes the in-trial/druggable cell from 5.97e-05 to 2.39e-04, dropping it a significance tier.

For the specificity test, Holm correction is applied separately within each panel rather than
pooled across both panels' six comparisons, since each panel's three comparisons share one
positive set; pooling all six under one Holm correction would demote two of the six significance
results by one tier.

For Figure 5, BH is applied as a single family across all eight enrichment tests it displays:
panel b's four depth-50 values (two arms by two stage thresholds) and panel c's four off-depth
phase-I values. Panel b and panel c share the depth-50 bar, so correcting each panel separately
would assign that bar two different p-values. BH is used rather than Holm because the analysis is
an FDR-appropriate enrichment screen, and at this family size Holm would push both genetics-arm
values above 0.05. Phase-III columns of the off-depth sweep are written but never plotted and are
excluded from the family, so they carry no `_bh` column.

For the leakage-controlled arms, the controlled and standard designs are treated as one
pre-specified family, since both report the same contrast under two evaluation designs on
Figure S3, and each arm's significance marker reflects the Holm-corrected p-value; the
published-value assertions checked in code remain pinned to the raw, uncorrected DeLong p-values.

### Other significance tests

- **DeLong test.** `build_leakage_controlled.py` computes a paired DeLong test on the single
  subframe where both methods' scores are defined under their own missing-value policies, so the
  two score vectors cannot be misaligned by independently filtering each method's gene set. In the
  standard arm the apparent GPS lead is not significant (p = 0.73).
- **Fisher exact test with Woolf 95% CI.** `make_enrichment_table.py` computes, for each of 4
  group x background cells, the odds ratio of IEI membership via scipy's `fisher_exact`, which
  returns the sample odds ratio (a*d)/(b*c), the estimator the Woolf standard error applies to,
  rather than a conditional MLE.
- **Mann-Whitney.** `build_score_by_group.py` runs 3 pairwise drug-status tests on the
  cross-fitted `pu_score_cv`, not the published `pu_score`. The approved genes are training
  positives on every PU-bagging iteration, so their published scores are in-sample; stage 2b of
  `pu_target_model.py` re-scores every gene with 5 folds over the positives at T=200, each positive
  from the one fold that held it out, so no comparison now rests on in-sample positive labels.
  One asymmetry survives cross-fitting, because it is not a positive-labelling artefact: `U_idx`
  is `role == 'unlabeled'`, so non-target genes are drawn as pseudo-negatives and in-trial genes
  never are. Out-of-bag averaging drops a gene's own negative iterations but not the boundary
  those iterations shaped, so the in-trial versus non-target contrast still sets a group that
  never entered training against the group that supplied every pseudo-negative.
  `build_supp_stat_sheets.py` checks that contrast by asserting the two U statistics sum to
  `n_it x n_nt`.
- **Bootstrap.** `build_immune_stages_data.py` draws each arm's control from that model's own
  ranking with its own top-200 excluded, using its own re-seeded `RandomState`; fold and p-value
  are computed against that arm's own null. `build_supp_stat_sheets.py` reports rates from 10,000
  bootstrap draws (`mean()*100` over 50 control genes per draw for panel b), with `n_ctrl_per_draw`
  fixed at 50 for panel b and varying with depth for panel c. `p_bootstrap` is the raw proportion
  of draws at or above the observed rate, not the continuity-corrected `(k+1)/(n+1)` estimator used
  in Tables 7 and 8. Recovery@k confidence intervals (`recovery_eval.py`) use a separate gene-level
  bootstrap: pool genes are resampled with replacement and hits@k recomputed at each resample to
  produce a percentile CI.
- **Permutation null.** `recovery_eval.py`'s null holds the ranking fixed and permutes which pool
  genes are labeled emergent, preserving the emergent count, then recomputes hits@k for an
  empirical p-value and fold-enrichment; this differs from a null that retrains on shuffled labels.
  `finalize_label_permutation.py` recomputes p-values against the current AUC ladder from an
  existing null table (`outputs/permutation/label_permutation_null.checkpoint.json` when present,
  else the committed `figure_data/label_permutation_null.csv`) without refitting the 1,000-draw
  null; all four `emp_p` values re-derive from the committed draws as `(k+1)/(n+1)`. The sheet
  reports levels only: the increment rows it once carried tested no rung against another, since
  label permutation drives every rung to chance, and the scrambled-feature control is that test. The
  discordance STRING test (Makefile) draws 500 matched nulls, giving an empirical p floor of 0.002;
  a reported p near that value reflects the floor rather than the true effect size, and lowering
  the draw count raises the floor and weakens the claim.
- **Class prior estimation.** `pu_target_model.py` estimates the Elkan-Noto class prior c as
  `E[score | positive]` using cross-validated held-out positive scores rather than in-bag
  optimistic scores; the estimated positive fraction in U is `E[score | U] / c`.

### Temporal holdout permutation resolution

`build_temporal_holdout.py`'s freeze-year sweep retrains both feature sets at each freeze year, so
the reported increment receives the same paired DeLong test as the headline result, using
n_perm=5000 permutations; the headline result itself uses n_perm=10000. The p-value for the same
freeze year therefore appears at two different permutation resolutions in the paper, both valid
floors on the true p-value. T0=2014 gives 134 emergent-test genes versus 41 at T0=2015; using
`known | dropped` rather than `known` alone as holdout genes shifts the T0=2015 significance test
across the p<0.05 threshold (0.0436 to 0.0507) at n_emergent=41, a sensitivity of the threshold
near a borderline sample size rather than a difference in design. The effect size is positive at
all five freeze years (+0.027 to +0.050); this, not a count of years reaching p<0.05, is the
finding supported across choices. T0=2018 is excluded from the sweep because too few post-freeze
years remain for emergent genes to accumulate, and the comparison is not significant
(DeLong p=0.074).

---

## Data sources and frozen inputs

### External data pins

| Source | Pin | Consequence of drift |
|---|---|---|
| STRING | v12.0 | Edge table built once from the bulk STRING flat files and cached as a symbol-pair CSV, so matched null draws are set lookups rather than one API call per draw. |
| GOA GAF, UniProt annotation-depth stream | unpinned, latest | Fetching on a different date changes the GO and dark-proteome numbers; exact reproduction requires the archived copies in `data/raw/discordance/` rather than a re-fetch. |
| Open Targets | release 26.06 | Association scores move between releases. The release is recorded in `METHODS.md`, which is not in this repository, rather than stamped into the committed drug data files, and is checked at runtime by `fetch_vignette_panelC_gwas.py --expect-release`, which aborts under any other release. |

`fetch_vignette_panelC_gwas.py` derives `gwas_score` as the `gwas_credible_sets` datasource score
for STAT4 x MONDO:0005046, aggregating evidence from 21 diseases carrying 115 credible sets at
release 26.06. The contributing set is fetched directly from the evidence endpoint with
`enableIndirect: true` rather than derived from each disease's ancestors list, since the two
values could diverge across releases; `ASSERT_TAGS` enforces that they agree at 26.06 and aborts
if a future release breaks that agreement. The parent-term score is separately checked against
the model matrix: the API returns 0.8311524998689545 and `data/perturbseq/pu/pu_model_matrix.parquet`
stores 0.831152 (six decimals), asserted equal within `FEATURE_TOL`. An earlier version of this
panel labeled the displayed rows as independent evidence corroborating the model; five of the six
are in fact carried by `gwas_credible_sets`, the same datasource behind the model's `gwas_score`
feature, and the sixth (pansclerotic morphea) by rare-variant datasources (ClinVar/eva, Genomics
England, UniProt variants) that map onto the model's IEI feature (STAT4 is listed in
`data/reference/IEI_gene_list.csv` with IEI = 1). The panel therefore reports the composition of
the MONDO:0005046 rollup rather than an independent check.

### Perturb-seq signature matrices

The knockdown-signature matrices under `data/perturbseq/knn/` are float32 memmaps rebuilt from
the S3-hosted perturb-seq `.h5ad` file, about 1.3 GB each, and are never committed.
`signatures_Stim48hr.npy` alone is a 7156 x 10282 matrix, 294 MB; committing it would consume a
large share of the Git LFS free bandwidth budget of 1 GB per month and throttle access after a
few clones. It is regenerated by `make knn-signatures` (`src/analysis/vignette/build_knn_signatures.py`),
which streams it from the public source in about 3 to 4 minutes. `build_vignette_panelD_knn.py`
is the only step that reads it; when it is absent the script skips the cosine computation with a
message and leaves `figure_data/vignette_panelD_knn.csv` and `vignette_knn_stats.json` at their
the committed values, so a clean clone still completes `make tables` and passes
`make verify-tables` without it. The row and column gene-name lists beside the matrices are
committed (110 KB total), so their identity can be checked without the matrices themselves.
Values reproduce exactly from the committed matrix: cosines of CD3E 0.355, LCK 0.328, CD3G 0.322,
CD2 0.306, PPP3CA 0.290, PRKCQ 0.283, and `knn_p95` = 0.234. The approved-drug and lead-status
columns in this table are hand-curated and not derived from any other source in the repository.

### GPS per-gene scores

`data/specificity/gps_per_gene.csv` and `data/specificity/gps_drug_labeled_genes.csv` are derived
from Supplementary Table 16 of the GPS paper by `build_gps_per_gene.py` (`make gps-inputs`).
Reproduction status differs by column:

| Column | Definition | Reproduces frozen file |
|---|---|---|
| `gps_max_overall` | max GPS across all 399 phecodes in the source table | exactly, for all 14,677 genes present in both files; the 4,825 genes absent from the frozen file are exactly the genes absent from the source table |
| `gps_drug_labeled_genes.csv` | genes with a non-null Open Targets Indication in the source table | exactly, 453 genes; Open Targets Clinical Phase selects the same 453 genes, SIDER Indication a different set of 287 |
| `gps_max_immune` | max GPS restricted to the 16 phecodes in `data/reference/immune_phecodes_curated.csv` | no: of 2,314 comparable genes, 2,221 match and 93 differ (for example ADAR, 0.880 frozen versus 1.286 here); 238 genes are scored only in the frozen file and 225 only in this derivation. The phecode set underlying the frozen column cannot be recovered from the current 16-phecode curated set. |

These files feed `build_leakage_controlled.py`, which produces the GPS versus PU comparison in
section 4.5a (GPS 0.641 vs PU 0.636 on the standard task; +0.053 for PU on the leakage-controlled
arm, paired DeLong z = 3.62, p = 2.9e-4). The leakage-controlled arm excludes the 453 genes used
to train GPS, leaving 535 positives that neither method was trained on. Reproducing the frozen
gene counts (19,162 for GPS, 18,878 for the leakage-controlled arm, 19,087 for missense-z)
requires filling missing GPS scores with 0, since GPS coverage is limited, while missing
missense-z rows are dropped, since absence there means the value is unmeasured. The 453-gene list
is shipped as a Table 11 sheet via `shutil.copyfile`, byte-for-byte, rather than through pandas,
so the shipped file's bytes match exactly the file the reported scores were computed from.

`panelA_auc_ladder.csv` is the only input to `build_coalition_tables.py` sourced from outside
`figure_data/`. It is committed via an explicit negation rule in `.gitignore`, so a clean clone
includes it, and is written by `attribution_decomposition.py` alongside that script's other,
gitignored scratch outputs.

### Drug status table

The `drugs` and `immune_indications` columns in the drug status table are sourced from the
per-drug evidence file rather than the per-gene rollup. The rollup caps these columns at 15 drugs
and 12 indications per gene: 44 genes are truncated under the drug cap (for example HRH1, with
`n_drugs` = 74 but only 15 drugs listed) and 290 genes exceed the indication cap. The rollup's
`n_drugs` and `furthest_stage` values are untruncated and are used to cross-check the table built
from the evidence file. `latest_trial_year` is not part of the output: `drug_trial_dates.json`
records only first-trial years per drug, and no committed data source records a true last-trial
date; that would require a live ClinicalTrials.gov query rather than the committed inputs this
script reads.

### Vignette metadata

`gwas_total` (183 GWAS Catalog hits, in `data/case_vignette/vignette_meta.json`) is a
hand-maintained, per-gene GWAS Catalog hit count carried over from a read-only seed file and
cannot be re-derived offline. It is distinct from Open Targets' GWAS disease rollup, computed as
`n_diseases_with_gwas_evidence` in `vignette_panelC_gwas_meta.json` (46 for STAT4); the two counts
must not be confused.

### Temporal holdout (T0) labels

`retrain_at_T0.py` defines the T0 retrain labels as follows. A gene enters the positive set P only
if it has an approval for an immune indication on or before T0; a gene in trials but not approved
by T0 is not a positive and instead belongs to `holdout_genes`. EMERGENT is the first
any-condition trial entry after T0, restricted to drugs not yet in clinic by T0; the trial can be
for any condition, not only immune indications, and the same T0 date defines both sides of the
split, so no separately chosen threshold enters the label. For the T0 <= 2014 vintage, 78
ribosomal protein genes are dated only via the MT-3724 trial join and enter the test set on that
basis. `holdout_genes` (genes in clinic by T0 without an approval) are excluded from the
pseudo-negative draw but are not part of the emergent test set. Only the label (trial-entry year)
is timestamped at T0; constraint, GWAS, and perturb-seq features use current data, which closes
label leakage but leaves feature leakage open.

### Perturbational functional genomics gene universe

The gene partition for the perturbational functional genomics figure is drawn from `pu_labels.csv`'s
`pu_role` field, independent of the drug tables used for Figure 2. Panels a, c, and d report
counts of distinct genes; panel b reports a count of rows. Each panel drops missing values
independently, so panels can cover different gene subsets. The `neglog10_adjp_Th*` columns are
display-only, a monotone transform of `|zscore_*|` used for the figure and excluded from the
model matrix. All four panels are restricted to the 19,502-gene universe defined by
`pu_labels.csv`.

### Feature dictionary provenance

Feature importances in `feature_dictionary.csv` are not recomputed by
`build_feature_dictionary.py`; they are taken from `pu_target_model.py`'s
`clf.feature_importances_`, accumulated over `T_BAG` bags, with mean and standard deviation
written to `outputs/model/pu_feature_importance_stability.csv`. Refitting a second XGBClassifier
to recompute them would reintroduce the reproducibility risk the xgboost version pin and the
`colsample_bytree` column-order rule are intended to prevent. Definitions in the same file are
sourced partly from the nine short labels in `build_discordance.py` and the six perturbational sub-block
names in `recompute_block_signal.py`; the remainder is drawn from `METHODS.md`, section 2, which is not in this repository.

### Files with no producer

A supplementary source with no producer script is never assigned a table number in the
supplementary package. `gps_drug_labeled_genes.csv` is frozen in this sense but reproducible in
the sense that matters: `build_gps_per_gene.py` (`make gps-inputs`) re-derives it byte-identically.

---

## Figure and table numbering

The `final-outputs` target assembles only the display items used in the manuscript, numbered as
the manuscript numbers them, splitting main and supplementary items out of the mixed, unorganized
contents of `final_plots/`.

### Supplementary tables

`supplementary_data/` is not assembled. Its former contents corresponded to published copies of
`figure_data/` tables that reached no numbered display item; those statistics now appear in
Supplementary Tables 7, 9 and 11, sourced directly from `figure_data/`.

Numbered supplementary tables are generated from `tools/supplementary_tables.tsv` by
`tools/package_supplementary_tables.py`, rather than from `cp` lines in the Makefile. A numbered
name is a manuscript fact that changes whenever a table is inserted or dropped, and the source
files cannot carry that fact: their names in `figure_data/` form the stage-1 contract read by
producers, `figure_data/PROVENANCE.md`, and `make verify-tables`. Renumbering therefore edits one
column in the manifest and nothing else. Filenames in `final_plots/supplementary/` cannot be used
for this tracking either, since that directory holds both build outputs and frozen artifacts used
by other lists.

`package_supplementary_tables.py` validates the manifest numbering for gaps and duplicates,
refuses the quarantined frozen GPS tables, and derives each output's column dictionary from that
file's actual header, never from the manifest description; a manifest description that does not
match the file is exposed by this check rather than trusted. It runs before
`write_submission_readme.py`, which lists the directory it creates.

Numbered supplementary tables are packaged directly from `figure_data/` and are never copied into
`final_plots/`, which keeps them out of `final_plots.sha256` and lets a new table be added without
a baseline refresh.

### Known panel/filename mismatches

- Figure 4 panel f (perturb-seq feature block signal, split by weak vs strong genetic prior) is
  computed by `recompute_block_signal.py` and written to
  `figure_data/panelE_block_signal.csv`; the file retains the `panelE` name though the manuscript
  panel is labeled f.
- Supplementary Table 15 (`add_knn_approved_drugs.py`) depends on signature matrices produced by
  `build_knn_nearest_target.py`; when those are absent both scripts produce no output and leave
  `outputs/knn/` empty, since the annotated tables are already committed in `figure_data/`. Its
  "approved" label reflects only `furthest_stage_drugs` rows with `furthest_stage == 'Approved'`,
  not the broader `drugs` column, which also includes assets still in trial.

---

## Reproduction boundary

### What `make submission` re-derives, and what it does not

`make submission` takes a fresh clone to `final_outputs/` in one command, offline, under the
pinned environments. It runs `check-versions`, `tables`, `figures` and `final-outputs` in that
order. What it re-derives, precisely:

| | Count | Status under `make submission` |
|---|---|---|
| `figure_data/` tables rewritten by `make tables` | 85 of 96 | Re-derived from committed inputs, then byte-compared by `verify-tables` |
| `figure_data/` tables belonging to opt-in targets | 11 of 96 | **Not re-derived.** Read at their committed values; `verify-tables` compares them against themselves and passes trivially |
| Rendered figures | 10 | Re-rendered from `figure_data/`; byte-compared by `verify-figures` and reported rather than enforced |
| Numbered supplementary tables | 15 | Packaged directly from `figure_data/` |

The 11 tables outside `make tables` are the 6 `discordance_*` caches, the 3
`label_permutation_*` files, and the 2 `vignette_panelC_gwas_*` files. Each has a producer in
this repository, but that producer needs the network or hours of compute and therefore sits
behind an opt-in target: `discordance-network`, `permutation-null`, `vignette-panelc`. The
`label_permutation_*` files were regenerated under the current pins on 2026-09-06 (see
**Permutation null** above). They remain the weakest point only in that no second run has
confirmed they re-derive.

Three further inputs are outside the offline path:

- **The perturb-seq signature matrices** (`data/perturbseq/knn/signatures_*.npy`, 281-294 MB
  each) are not committed. `make knn-signatures` streams them. Absent them, three stage-1
  steps skip and leave five committed `figure_data/` tables at their committed values;
  `make tables` prints exactly which ones at the end of the run, via
  `src/analysis/common/_skip_ledger.py`.
- **`data/raw/`** is not committed, and the fetch and subset scripts that would populate it
  are outside this reproduction subset. Every input the offline build needs is committed
  under `data/`, so this affects only the opt-in targets that read `data/raw/`.
- **Producerless files.** A small number of committed inputs have no producing script in this
  repository: the external lists in `data/prospective/`, two files in `data/specificity/`,
  and the hand-curated phecode list in `data/reference/immune_phecodes_curated.csv`.
  `figure_data/PROVENANCE.md` records which, and says `unknown` rather than guessing where
  the origin is not recorded.

Stated plainly: `make submission` reproduces every number the manuscript reports from the
committed inputs, and it does not re-derive those committed inputs from their primary sources.

### Three levels of reproduction

The repository draws a line between three levels of reproduction, and the checks below cover different sides of it: the fresh-clone stage-2 path (`verify-tables`, `verify-figures`), the full re-derivation path (`make tables`, `make figures`, `make final-outputs`), and the opt-in producers that regenerate committed inputs from external sources (`permutation-null`, `discordance-network`, `vignette-panelc`, `gps-inputs`, `knn-signatures`).

### What "reproduced" means for `verify-tables`

`figure_data/` is committed, so git already records the reference state; a clean working tree after a rebuild is the reproduction check itself, and no separate baseline copy of the tables is kept. `tools/write_baseline.sh` classifies files the same way, by asking git rather than maintaining an independent file list. The check requires a git checkout and does not work from a source tarball. The gate excludes `.md` files from the diff itself, not only from the reported count, which is computed with `grep -cv '\.md$'`; this keeps edits to documentation such as `PROVENANCE.md`, which lives in the same directory as the tables, from being reported as a table reproduction failure.

### Baseline generation

`tools/write_baseline.sh <py> <rscript>` produces `final_plots.sha256` from a clean clone built with the pinned interpreters, only after `make verify` has been reasoned about. It defines the reference against which every later verification is judged, so a baseline built from a stale or mismatched environment becomes silently authoritative for all subsequent checks.

- The hashing walk is recursive and null-delimited. A non-recursive listing such as `shasum -a 256 final_plots/*` omits nested outputs under `final_plots/supplementary/`; a space-delimited walk also splits filenames containing spaces.
- The Makefile exports the R library path for every R step. Run standalone, as documented, an ordinary shell resolves R packages from the system library instead of `.rlib/`; since the baseline defines the reference for all later verification, stamping the wrong package version causes `check_r_versions.R` to mark an otherwise correct tree as PROVISIONAL. The script resolves packages from `.rlib/` to match the environment that rendered the hashed figures.
- `git rev-parse HEAD` names a commit while the hashes come from the working tree; the script confirms the two agree before stamping a baseline against uncommitted work.
- Only build-relevant paths are covered. `figure_data/` is included even though `make tables` leaves it dirty, since stage 2 reads it and uncommitted tables there change the rendered bytes. `drafts/`, `logs/`, and `final_plots/` itself are excluded: the deliverable tree is the thing being hashed and is expected to be dirty, and the other two feed nothing into the build.

### Build scratch and cleanup scope

| Path | Committed | Removed by `make clean` | Notes |
|---|---|---|---|
| `outputs/` | No (gitignored) | n/a, always scratch | Everything stage 1 emits that no figure reads. A fresh clone reruns `make tables` to obtain gene and PU scores, since none are cached. Each negation in the `.gitignore` rule corresponds to a file a retained stage-1 script still writes; removing one leaves an untracked file in `git status` after `make tables`. |
| `data/raw/discordance/` | No | No | About 135 MB of one-time external downloads, reused only by `make discordance-network`; re-fetching costs roughly 10 minutes of network time. |
| `final_outputs/` | No | Yes | Rebuilt by `make final-outputs` as a byte-identical copy of `final_plots/` plus two generated READMEs (the top-level README and the numbered tables' column dictionary), performing no computation. |
| `figure_data/` (numbered supplementary tables) | Yes | No | Needs no separate clean entry: it is never copied into `final_plots/`. |

### Cross-figure consistency (Panel D recovery)

Panel D's boundary of reproduction extends beyond its own script to two other producers, and agreement is by construction rather than by re-check:

- The immune-indication membership rule (furthest immune-indication stage at or above Phase 1) matches `build_immune_stages_data.py`, so Panel D and the trial-validation figure agree by construction.
- The full-model bar uses the production PU ranking, `outputs/model/full_model_pu_scores.csv`, identical to Figure 5b; its top-50 count matches Figure 5b's count.
- The genetics-only bar is re-scored inside `make_panelD_recovery.py` because the production run does not persist a genetics-only ranking. `SEED` and `T_BAG` are imported from `pu_target_model` rather than restated, so the bar tracks the production configuration.
- Results are bit-reproducible under the pinned xgboost version.

Genetics-only ranking persistence, in the same script, writes all 19,502 rows rather than only the ranked pool, so the pool filter (`pu_role != 'P'`) can be audited rather than assumed. The rank column is stored explicitly because `pandas.sort_values` defaults to an unstable quicksort, so tie order among equal scores is not reproducible unless the resulting rank is persisted.

---

## Conventions

- Python scripts locate the repository root by walking up from `__file__`, never from the current working directory.
- Comments carry the reasoning behind a decision - why an order is load-bearing, what a previous version got wrong. That density is the house style across the analysis and figure scripts; it is preserved here rather than stripped.
- Counts are computed rather than hardcoded where possible; hardcoded integers in this document have drifted before. The same principle governs naming: where a string serves as both a data key and a printed label - for example the group-name strings used as keys in `score_by_group_stats.json` (written by `build_score_by_group.py`) and as the labels printed in the Supp Table 7 sheets (`build_supp_stat_sheets.py`) - the names are read through a single shared list rather than hardcoded at each use site, so a rename of the JSON keys without a matching sheet update raises a KeyError instead of silently producing a sheet of NaNs.
- Typography: the figure set is unbolded. No `face = "bold"` appears in any renderer, panel tags included; mixed weight is visible in the assembled PDF.
- Build outputs are produced by the Makefile and are not committed; a clean clone regenerates them. `final_plots/` mixes frozen artifacts that have no producing script with generated files at the same directory level, so its `.gitignore` rules are written per file rather than per directory. The `.gitkeep` placeholders inside ignored directories are tracked with `git add -f`, since a `!.gitkeep` negation cannot re-include a file whose parent directory is excluded.
- `make tables` starts with `rm -rf logs`; per-step logs are written to `logs/NN_<script>.log`, with a summary in `logs/_summary.tsv`.
- A dirty `git status` after `make tables` is expected, since the step rewrites the committed contract in place; `make verify` is the signal to check, not `git status` on its own.

### Running one step

There is no test runner: the unit of work is one script. To rerun a single step exactly as the build does:

```bash
PY=.venv/bin/python RSCRIPT=Rscript tools/run_step.sh py src/analysis/discordance/build_discordance.py
PY=.venv/bin/python RSCRIPT=Rscript tools/run_step.sh R  src/figures/fig_discordance.R
```

Or directly (the Makefile exports `PYTHONPATH=src`):

```bash
PYTHONPATH=src .venv/bin/python src/analysis/specificity/build_specificity_matrix.py
R_LIBS_USER=$PWD/.rlib Rscript src/figures/fig_specificity.R
