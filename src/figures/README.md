# src/figures — stage 2

10 R renderers. Each reads `figure_data/` and nothing else, and writes its deliverable directly into `final_plots/`. Every count, star, p-value and threshold is precomputed by stage 1, so a figure needing a new number needs a stage-1 change.

Run with `make figures` (~15 s, no model fitting).

# Renderers

| Script | Manuscript figure | Output under `final_plots/` |
|---|---|---|
| `fig_genetic_data.R` | Figure 2, genetic priors | `genetic_data.png` |
| `fig_perturbational_functional_genomics.R` | Figure 3, perturbational functional genomics | `perturbational_functional_genomics_with_thde.png` |
| `fig_model_analysis.R` | Figure 4, model analysis composite | `figure_model_analysis_composite.png` |
| `fig_immune_stages.R` | Figure 5, trial validation | `figure_trial_validation_immune_stages_R.png` |
| `fig_temporal_holdout.R` | Figure 6, prospective temporal holdout | `figure_prospective_temporal_holdout.png` |
| `fig_specificity.R` | Figure 7, specificity control | `specificity_immune_vs_cardiac_faceted.png` |
| `fig_stat4_vignette.R` | Figure 8, STAT4 vignette | `figure_stat4_vignette.png` |
| `fig_label_permutation.R` | Supplementary Figure 1 | `figure_label_permutation.png` |
| `fig_discordance.R` | Supplementary Figure 2 | `supplementary/figure_discordance_genetics_vs_full.png` |
| `fig_leakage_controlled.R` | Supplementary Figure 3 | `supplementary/leakage_controlled_comparison.png` |

- `tools/supplementary_figures.tsv` is the single source for the S-numbers.
- `make supplementary` re-renders only the three supplementary renderers, listed as `SUPP_RENDERERS_R` in the Makefile, and writes to its own `supplementary/` tree.
- Each renderer writes a `.png`. On macOS it also writes a vector `.pdf` to `final_plots/pdf/`; other platforms produce the PNG only.
- Numbered supplementary tables are packaged from `figure_data/` by `tools/package_supplementary_tables.py` and never pass through `final_plots/`.

# Shared helpers

| File | Provides |
|---|---|
| `palette.R` | `TARGET_COLORS`, the drug-status colour mapping shared across figures |
| `vector_output.R` | `save_figure(plot, out_png, out_pdf, width_cm, height_cm)` — PNG through the ragg AGG device, PDF through quartz (macOS only) |

# Environment

- R 4.4.1 and the package pins in `R-requirements.txt`, checked by `tools/check_r_versions.R` (`make check-versions`).
- `make setup` provisions them into `.rlib/`, which the Makefile adds to the R library path when present.
- `systemfonts` and `textshaping` determine glyph layout for the ragg device and are pinned for that reason. `ggtext` and `gridtext` back `fig_stat4_vignette.R` panel c.
- Exact versions and system-library notes are in `REPRODUCIBILITY.md`.

# Verification

`make verify-figures` compares `final_plots/` against `final_plots.sha256`. Differing PNGs and PDFs are reported but do not fail the gate; it fails when a baseline file is missing or when `final_plots/` contains a file the baseline does not list. Under the pinned R stack PNG bytes reproduce on the same machine; a differing PNG usually means the render did not use `.rlib/`. PDFs carry a creation timestamp, so they differ on every render.
