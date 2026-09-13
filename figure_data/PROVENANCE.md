# figure_data/ — provenance

Contract between the two build stages: `make tables` writes this directory, and `make figures` reads only from it. It is committed to git so that `make figures` runs from a fresh clone with no model fitting.

- One line per file: the stage-1 script that writes it.
- Scripts live under `src/analysis/`; see `src/analysis/README.md` for the folder map and `REPRODUCIBILITY.md` for the step order.
- 98 committed non-`.md` files: 87 written by `make tables`, 11 by opt-in targets. Verify the count with `git ls-files figure_data | grep -vc '\.md$'`.
- Opt-in files are written by targets outside `make tables`, so `make verify-tables` confirms they are unchanged rather than re-derived.
- `make tables` rewrites the stage-1 files in place, leaving the working tree dirty here; `make verify-tables` checks that they re-derive.

# Written by stage 1 (87 files)

## Genetic priors — Fig 2, `final_plots/genetic_data.png`

| File | Written by |
|---|---|
| `iei_enrichment_forest.csv` | `make_enrichment_table.py` |
| `iei_enrichment_axis.csv` | `make_enrichment_table.py` |
| `gwas_by_drug_status.csv` | `build_genetic_data_tables.py` |
| `mis_z_by_drug_status.csv` | `build_genetic_data_tables.py` |
| `genetic_data_group_n.csv` | `build_genetic_data_tables.py` |
| `gwas_zero_inflation.csv` | `build_genetic_data_tables.py` |

## Perturbational functional genomics — Fig 3, `final_plots/perturbational_functional_genomics_with_thde.png`

| File | Written by |
|---|---|
| `th_de_long.csv` | `build_perturbational_fg_tables.py` |
| `th_de_group_n.csv` | `build_perturbational_fg_tables.py` |
| `regulator_residual.csv` | `build_perturbational_fg_tables.py` |
| `regulator_residual_group_n.csv` | `build_perturbational_fg_tables.py` |
| `cytokine_counts.csv` | `build_perturbational_fg_tables.py` |
| `cytokine_group_n.csv` | `build_perturbational_fg_tables.py` |
| `cytokine_receptor_counts.csv` | `build_perturbational_fg_tables.py` |
| `cytokine_receptor_group_n.csv` | `build_perturbational_fg_tables.py` |

## Model analysis — Fig 4, `final_plots/figure_model_analysis_composite.png`

| File | Written by |
|---|---|
| `panelA_auc_ladder.csv` | `reshape_attribution_to_panels.py` |
| `panelA_significance.csv` | `reshape_attribution_to_panels.py`, rewritten in place by `compute_nb_significance.py` |
| `panelB_shapley_shares.csv` | `reshape_attribution_to_panels.py` |
| `panelC_marginal.csv` | `reshape_attribution_to_panels.py`, rewritten in place by `compute_nb_significance.py` |
| `panels_ABC_meta.json` | `reshape_attribution_to_panels.py` |
| `panelD_recovery.csv` | `make_panelD_recovery.py` |
| `panelE_block_signal.csv` | `recompute_block_signal.py` |
| `score_by_group.csv` | `build_score_by_group.py` |
| `score_by_group_stats.json` | `build_score_by_group.py` |

## Trial validation — Fig 5, `final_plots/figure_trial_validation_immune_stages_R.png`

| File | Written by |
|---|---|
| `immune_stages_top50.csv` | `build_immune_stages_data.py` |
| `immune_stages_counts.json` | `build_immune_stages_data.py` |

## Prospective temporal holdout — Fig 6, `final_plots/figure_prospective_temporal_holdout.png`

| File | Written by |
|---|---|
| `t0_split.json` | `build_temporal_holdout.py` |
| `scores_at_T0_2014.csv` | `build_temporal_holdout.py` |
| `auc_emergent.csv` | `build_temporal_holdout.py` |
| `t0_sensitivity.csv` | `build_temporal_holdout.py` |
| `t0_panel_stats.csv` | `build_temporal_holdout.py` |
| `emergent_targets_T0_2014.csv` | `build_temporal_holdout.py` |

## Specificity control — Fig 7, `final_plots/specificity_immune_vs_cardiac_faceted.png`

| File | Written by |
|---|---|
| `specificity_matrix.csv` | `build_specificity_matrix.py` |
| `specificity_panel_n.csv` | `build_specificity_matrix.py` |
| `delong_specificity_tests.csv` | `compute_delong_specificity.py` |

## STAT4 vignette — Fig 8, `final_plots/figure_stat4_vignette.png`

| File | Written by |
|---|---|
| `vignette_meta.json` | `build_vignette_meta.py` |
| `vignette_panelB_axes.csv` | `build_vignette_panelB.py` |
| `vignette_axes_supp.csv` | `build_vignette_panelB.py` |
| `vignette_score_hist.csv` | `pu_target_model.py` |
| `vignette_panelD_knn.csv` | `build_vignette_panelD_knn.py` |
| `vignette_knn_stats.json` | `build_vignette_panelD_knn.py` |

## Label permutation and scrambled-feature control — S1, `final_plots/figure_label_permutation.png`

| File | Written by |
|---|---|
| `scrambled_control_stats.csv` | `build_supp_stat_sheets.py` |

## Genetics-vs-full discordance — S2, `final_plots/supplementary/figure_discordance_genetics_vs_full.png`

| File | Written by |
|---|---|
| `discordance_scatter.csv` | `build_discordance.py` |
| `discordance_drivers.csv` | `build_discordance.py` |
| `discordance_features.csv` | `build_discordance.py` |
| `discordance_validation.csv` | `build_discordance.py` |
| `discordance_core_genes.csv` | `build_discordance.py` |
| `discordance_stats.json` | `build_discordance.py` |

## Leakage-controlled comparison — S3, `final_plots/supplementary/leakage_controlled_comparison.png`

| File | Written by |
|---|---|
| `leakage_controlled_comparison.csv` | `build_leakage_controlled.py` |
| `leakage_controlled_delong.json` | `build_leakage_controlled.py` |
| `gps_drug_labeled_genes.csv` | `build_leakage_controlled.py` |

## Ranked atlas — Supp Table 1

| File | Written by |
|---|---|
| `ranked_atlas.csv` | `build_ranked_atlas_table.py` |

## Feature dictionary — Supp Table 2

| File | Written by |
|---|---|
| `feature_dictionary.csv` | `build_feature_dictionary.py` |

## Drug-status annotations — Supp Table 3

| File | Written by |
|---|---|
| `drug_status_annotations.csv` | `build_drug_status_table.py` |

## Univariate distributions by drug status — Supp Table 4

| File | Written by |
|---|---|
| `univariate_iei_enrichment.csv` | `build_univariate_by_group.py` |
| `univariate_gwas.csv` | `build_univariate_by_group.py` |
| `univariate_gwas_pct_gt0.csv` | `build_univariate_by_group.py` |
| `univariate_gwas_tests.csv` | `build_univariate_by_group.py` |
| `univariate_mis_z.csv` | `build_univariate_by_group.py` |
| `univariate_mis_z_tests.csv` | `build_univariate_by_group.py` |
| `univariate_th_de.csv` | `build_univariate_by_group.py` |
| `univariate_th_de_fdr_enrichment.csv` | `build_univariate_by_group.py` |
| `univariate_residual.csv` | `build_univariate_by_group.py` |
| `univariate_residual_tests.csv` | `build_univariate_by_group.py` |
| `univariate_cytokines.csv` | `build_univariate_by_group.py` |
| `univariate_cytokine_receptors.csv` | `build_univariate_by_group.py` |
| `univariate_cytokine_distribution.csv` | `build_univariate_by_group.py` |
| `univariate_cytokine_regulator_fraction.csv` | `build_univariate_by_group.py` |
| `univariate_group_n.csv` | `build_univariate_by_group.py` |

## Orthogonality of functional genomics to the genetic priors — Supp Table 5

| File | Written by |
|---|---|
| `orthogonality_genetic_vs_fg.csv` | `build_orthogonality_tables.py` |
| `orthogonality_genetic_pairs.csv` | `build_orthogonality_tables.py` |

## Evidence coalitions and attribution — Supp Table 6

| File | Written by |
|---|---|
| `coalition_auc.csv` | `build_coalition_tables.py` |
| `coalition_ladder.csv` | `build_coalition_tables.py` |
| `coalition_significance.csv` | `build_coalition_tables.py` |
| `coalition_attribution.csv` | `build_coalition_tables.py` |
| `coalition_strata.csv` | `build_coalition_tables.py` |

## Held-out in-trial recovery — Supp Table 7

| File | Written by |
|---|---|
| `heldout_trial_roc_curves.csv` | `build_heldout_trial_auc.py` |
| `heldout_trial_auc.json` | `build_heldout_trial_auc.py` |
| `heldout_trial_stats.csv` | `build_supp_stat_sheets.py` |
| `score_by_group_summary.csv` | `build_supp_stat_sheets.py` |
| `score_by_group_tests.csv` | `build_supp_stat_sheets.py` |

## Flattened stats sheets — Supp Tables 10, 11 and 14

| File | Written by |
|---|---|
| `immune_stages_bootstrap.csv` | `build_supp_stat_sheets.py` |
| `immune_stages_depth_sweep.csv` | `build_supp_stat_sheets.py` |
| `leakage_controlled_delong_stats.csv` | `build_supp_stat_sheets.py` |
| `vignette_knn_summary.csv` | `build_supp_stat_sheets.py` |

## Nearest approved-target knockdown neighbour — Supp Table 15

| File | Written by |
|---|---|
| `knn_nearest_target_Rest.csv` | `add_knn_approved_drugs.py` |
| `knn_nearest_target_Stim8hr.csv` | `add_knn_approved_drugs.py` |
| `knn_nearest_target_Stim48hr.csv` | `add_knn_approved_drugs.py` |

# Written by opt-in targets (11 files)

## `make discordance-network` (6 files)

| File | Written by |
|---|---|
| `discordance_network.csv` | `fetch_discordance_network.py` |
| `discordance_network_edges.csv` | `fetch_discordance_network.py` |
| `discordance_network_nodes.csv` | `fetch_discordance_network.py` |
| `discordance_go_enrichment.csv` | `fetch_discordance_network.py` |
| `discordance_dark_proteome.csv` | `fetch_discordance_network.py` |
| `discordance_dark_proteome_axes.csv` | `fetch_discordance_network.py` |

## `make permutation-null` (3 files)

| File | Written by |
|---|---|
| `label_permutation_null.csv` | `run_label_permutation.py` |
| `label_permutation_pvalues.csv` | `run_label_permutation.py`, or `finalize_label_permutation.py` to recompute p-values alone |
| `scrambled_feature_control.json` | `scrambled_feature_control.py` |

## `make vignette-panelc` (2 files)

| File | Written by |
|---|---|
| `vignette_panelC_gwas_diseases.csv` | `fetch_vignette_panelC_gwas.py` |
| `vignette_panelC_gwas_meta.json` | `fetch_vignette_panelC_gwas.py` |

# Inputs read from outside this directory

| Path | What it is |
|---|---|
| `data/perturbseq/pu/pu_model_matrix.parquet` | the 19,502 x 28 feature matrix (24 features + 4 meta); the committed starting point for stage 1 |
| `data/perturbseq/pu/pu_labels.csv` | the label table and `pu_role` partition |
| `data/perturbseq/knn/signatures_*.npy` | knockdown signature matrices; not committed, streamed by `make knn-signatures`. The row and column gene-name lists beside them are committed |
| `data/specificity/gps_per_gene.csv`, `gps_drug_labeled_genes.csv` | GPS inputs to `build_leakage_controlled.py`; committed, re-derived by `make gps-inputs` |
| `data/reference/immune_phecodes_curated.csv` | the 16 curated immune phecodes used by `build_gps_per_gene.py`; hand-curated, no producing script |
| `data/reference/go_term_labels.json` | GO-id-to-term-name map; extended in place by `fetch_discordance_network.py` |
| `data/raw/discordance/` | STRING and GO archives read by `make discordance-network`; not committed |
| `outputs/model/full_model_pu_scores.csv` | stage-1 scratch written by `pu_target_model.py` and read by four later steps; not committed |
| `outputs/model/genetics_only_pu_scores.csv` | stage-1 scratch written by `make_panelD_recovery.py`; not committed |
| `outputs/model/pu_feature_importance_stability.csv` | bagging importances written by `pu_target_model.py`, read by `build_feature_dictionary.py`; not committed |
| `outputs/knn/` | stage-1 scratch written by `build_knn_nearest_target.py`, annotated by `add_knn_approved_drugs.py`; not committed |

See `DATA_AVAILABILITY.md` for the external provenance and redistribution terms of every committed input.
