# src/figures — stage 2 renderers

## Renderers

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

`tools/supplementary_figures.tsv` is the single source for the S-numbers. `make supplementary`
re-renders only the three supplementary renderers, listed as `SUPP_RENDERERS_R` in the Makefile.

Each renderer writes a `.png`. On macOS it also writes a vector `.pdf` to `final_plots/pdf/`;
other platforms produce the PNG only.

## Shared helpers

| File | Provides |
|---|---|
| `palette.R` | `TARGET_COLORS`, the drug-status colour mapping shared across figures |
| `vector_output.R` | `save_figure(plot, out_png, out_pdf, width_cm, height_cm)`, which writes the PNG through the ragg AGG device and the PDF through a quartz device (macOS only) |

## Environment

The R package pins are in `R-requirements.txt` and are checked by `tools/check_r_versions.R`
(`make check-versions`). `make setup` provisions them into `.rlib/`, which the Makefile adds to
the R library path when present. Font handling matters here: `systemfonts` and `textshaping`
determine glyph layout for the ragg device, and are pinned for that reason.

## Figure verification

`make verify-figures` compares `final_plots/` against `final_plots.sha256`. Under the pinned R
stack PNG bytes reproduce; a differing PNG usually means the render did not use `.rlib/`, or ran
on a machine with a different font and rasterization stack. PDFs carry a creation timestamp, so
they differ on every render. Differing PNGs and PDFs are reported but do not fail the gate. It
fails when a baseline file is missing or when `final_plots/` contains a file the baseline does not
list. See REPRODUCIBILITY.md, `make verify-tables` and `make verify-figures`: the five-class system.
