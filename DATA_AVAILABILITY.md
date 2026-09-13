# Data Availability

Provenance and redistribution terms for every external input. Per-file provenance for the derived tables is in `figure_data/PROVENANCE.md`.

- Most external data is committed under `data/` and is sufficient to run `make tables` and `make figures` offline from a clean clone.
- Three perturb-seq signature matrices and everything under `data/raw/` are not committed, owing to size, and are obtained by opt-in Makefile targets.
- The code is MIT-licensed. The redistributed data files are not, and remain under their original providers' terms.

# External sources

| Source | What is used | Where it lives | Access |
|---|---|---|---|
| **Open Targets**, release 26.06 (platform.opentargets.org) | approved/in-trial drug-by-gene and drug-by-evidence tables; immune-disorder GWAS and burden scores (MONDO:0005046, `gwas_score`, `gene_burden_score`); genetic-association scores for immune and cardiac indications | `data/data_drug/all_drugs_approved_and_in_trial_by_gene_ot.csv`, `data/data_drug/all_drugs_approved_and_in_trial_evidence_ot.csv`, `data/data_drug/approved_target_genes.txt`, `data/specificity/cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv`, `data/perturbseq/pu/pu_model_matrix.parquet` | committed. `make vignette-panelc` re-queries the live API and asserts release 26.06 |
| **Perturb-seq**, genome-scale CRISPRi in primary human CD4+ T cells; Zhu, Dann et al. (2025) | knockdown `zscore` signatures (log2FC / lfcSE); PU-model feature inputs; cytokine and cytokine-receptor significance counts | `data/perturbseq/knn/` (gene-name lists committed, matrices not), `data/perturbseq/pu/pu_model_matrix.parquet`, `data/perturbseq/pu/pu_labels.csv`, `data/perturbseq/sources/cytokine_significant_counts.csv`, `data/perturbseq/sources/cytokine_receptor_significant_counts.csv` | released via the CZI Virtual Cells Platform at `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`. Must be cited in any work using these signatures |
| **T-helper polarization DE** (github.com/emdann/GWT_perturbseq_analysis_2025) | DE z-scores of polarized T-helper subsets against resting Th0 | `data/perturbseq/pu/pu_model_matrix.parquet` (`zscore_Th1`, `zscore_Th2`, `zscore_Th17`, `zscore_Treg`) | committed as feature columns |
| **HGNC** complete gene set, accessed 2026-07-16 | `symbol`, `status`, `alias_symbol`, `prev_symbol` | `data/perturbseq/sources/hgnc_symbol_subset.txt` | committed as a 4-of-54 column subset; see below |
| **gnomAD** v4.1.0 | missense constraint (`mis.z_score`), LOEUF (`lof.oe_ci.upper`) | `data/perturbseq/pu/pu_model_matrix.parquet` | committed as feature columns |
| **IUIS** inborn errors of immunity list, 2024 | gene-level IEI status | `data/reference/IEI_gene_list.csv` | committed |
| **Druggable genome**, Finan et al. 2017 (doi:10.1126/scitranslmed.aag1166) | gene-level druggability annotation | `data/reference/druggable_genome_gene_list.csv` | committed |
| **GPS** (Genetic Priority Score), zenodo.org/records/10095684 | per-gene, per-phecode priority scores | `data/specificity/gps_per_gene.csv`, `data/specificity/gps_drug_labeled_genes.csv` committed; source workbook not | source fetched by `make gps-inputs` |
| **Immune phecode list** | 16 immune phecodes restricting GPS scores to immune indications | `data/reference/immune_phecodes_curated.csv` | committed; hand-curated, no producing script |
| **Drug trial dates** (ChEMBL, ClinicalTrials.gov) | first-trial year per drug, approval year per gene, for the temporal holdout labels | `data/prospective/drug_trial_dates.json`, `data/prospective/positives_dated.csv` | committed; no producing script |
| **GWAS Catalog / Open Targets Genetics** | per-gene GWAS Catalog hit count; Open Targets GWAS-evidence disease rollup | `data/case_vignette/vignette_meta.json`, `figure_data/vignette_panelC_gwas_meta.json` | the hit count is hand-maintained with no producing script; the rollup is re-derived by `make vignette-panelc` |
| **STRING** v12.0 (string-db.org) | protein-protein interaction edges (protein links, `protein.info`, `protein.aliases`) | not committed; fetched to `data/raw/discordance/` by `make discordance-network`. Derived tables committed under `figure_data/` | raw archives not redistributed |
| **Gene Ontology** (geneontology.org) | UniProt-GOA GAF and a UniProt annotation-depth table | not committed; fetched to `data/raw/discordance/` by `make discordance-network`. Derived tables committed under `figure_data/` | raw files not redistributed; release unpinned (latest at fetch time) |

Two further committed inputs are gene-universe tables assembled from the sources above rather than separate providers:

| File | Shape | Read by |
|---|---|---|
| `data/full_gene_list.tsv` | 19,502 x 17 | `make_enrichment_table.py` (the `all_genes` background), `build_genetic_data_tables.py` |
| `data/full_gene_list.with_perturbseq.tsv` | 19,502 x 37 | `build_perturbational_fg_tables.py`, `build_vignette_panelB.py`, `build_gps_per_gene.py` |

`data/reference/go_term_labels.json` is a committed GO-id-to-term-name map, extended in place by `fetch_discordance_network.py` from QuickGO for ids it cannot resolve locally.

# Not committed

## Perturb-seq signature matrices

Float32 memmaps rebuilt from the perturb-seq `.h5ad` source:

| File | Size | Dimensions |
|---|---|---|
| `data/perturbseq/knn/signatures_Rest.npy` | 281 MB | 6,838 x 10,282 |
| `data/perturbseq/knn/signatures_Stim8hr.npy` | 292 MB | 7,112 x 10,282 |
| `data/perturbseq/knn/signatures_Stim48hr.npy` | 294 MB | 7,156 x 10,282 |

```
make knn-signatures     # network, ~3-4 min
```

- Rows are perturbed genes with a significant on-target knockdown; columns are the 10,282 measured genes; values are the `zscore` layer (log2FC / lfcSE).
- Streams the `zscore` layer over HTTP range requests from `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`. Rows are de-duplicated and kept where `ontarget_significant` is true and `n_cells_target >= 50`.
- Read downstream by `build_vignette_panelD_knn.py` (`Stim48hr` only) and `build_knn_nearest_target.py` (all three).
- The accompanying row and column gene-name lists (`signatures_*.genes.txt`, `var_gene_names.txt`, 110 KB total) are committed, so row and column identity is checkable without the matrices.
- Not committed because Git LFS's 1 GB/month free bandwidth budget would cover about three clones of the largest file.

## data/raw/

Gitignored (`/data/raw/*`) and absent on a clean clone:

- `data/raw/discordance/` — STRING v12.0 protein links (~83 MB) plus `protein.info` and `protein.aliases`, the UniProt-GOA GAF (~11 MB), and a UniProt annotation-depth table (~3.5 MB); ~135 MB total. Obtained with `make discordance-network`; reused on later runs and not removed by `make clean`.
- `data/raw/gps/GPS_allgenes.xlsx` — the GPS source workbook (~5.6 MB). Obtained with `make gps-inputs`.

# The HGNC file

`data/perturbseq/sources/hgnc_symbol_subset.txt` is a column subset of the upstream HGNC complete set, not the verbatim release.

- Retains 4 of 54 columns (`symbol`, `status`, `alias_symbol`, `prev_symbol`) for all 45,021 rows — the only columns any consumer here reads. Reuse needing other columns requires the complete file from HGNC.
- HGNC does not version its complete set, so the release is identified by content: 45,021 approved records, latest `date_modified` 2026-07-14 (md5 `08386692953902f3880b52ba7d8b2f02` for the 54-column upstream file). `load_hgnc` prints both on every load.
- HGNC renames symbols between releases. `IL1RL1`, `IL1RL2` and `C17orf99` are approved in this release and become `IL33R`, `IL36R` and `IL40` later; `IL1RL2` is an approved drug target, so harmonizing against a different release re-keys it out of every join on `gene`.
- `load_hgnc` refuses to fetch the current release unless called with `allow_download=True`; a missing pinned file raises rather than falling back to the network.

# Files with no producing script

- the external lists in `data/prospective/`
- `data/specificity/ot_minus_clinical_per_gene.csv` and the GWAS Catalog hit count in `data/case_vignette/vignette_meta.json`
- `data/reference/immune_phecodes_curated.csv`

# Licensing

- The code is MIT-licensed; see `LICENSE`.
- Redistributed third-party data files remain under their original providers' terms; this repository does not relicense them. Anyone reusing or redistributing Open Targets, HGNC, STRING, Gene Ontology / UniProt-GOA, GWAS Catalog, the perturb-seq dataset of Zhu, Dann et al. (2025), gnomAD, the GPS dataset, or any other source named above should consult that provider's terms and cite the corresponding publication independently of this repository.
- Portions of the perturb-seq feature-engineering approach derive from emdann/GWT_perturbseq_analysis_2025 (MIT License, (c) 2025 Emma Dann). The implementing scripts are not part of this repository; only their committed output, `pu_model_matrix.parquet`, is.
