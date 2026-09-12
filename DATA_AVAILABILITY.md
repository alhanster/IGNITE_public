# Data Availability

This file is the single point of reference for the provenance and redistribution terms of every
external input used by this repository. Per-file provenance for the tables derived from those
inputs is recorded separately in `figure_data/PROVENANCE.md`; this file addresses the inputs
themselves.

Most external data is committed under `data/` and is sufficient to run `make tables` and
`make figures` offline from a clean clone. Three perturb-seq signature matrices and everything
under `data/raw/` are not committed, owing to size, and are obtained on demand by opt-in Makefile
targets described below. The code in this repository is MIT-licensed; the redistributed data
files are not, and remain under the terms of their original providers (see the licensing section
at the end of this file).

## External data sources

| Source | What is used | Where it lives in this repo | Access |
|---|---|---|---|
| Open Targets (release 26.06, https://platform.opentargets.org/) | Drug and association evidence: approved/in-trial drug-by-gene and drug-by-evidence tables, immune system disorder GWAS and burden scores (MONDO:0005046, `gwas_score`, `gene_burden_score`), and genetic-association scores for immune and cardiac indications | `data/data_drug/all_drugs_approved_and_in_trial_by_gene_ot.csv`, `data/data_drug/all_drugs_approved_and_in_trial_evidence_ot.csv`, `data/data_drug/approved_target_genes.txt`, `data/specificity/cardiac_all_drugs_approved_and_in_trial_by_gene_ot.csv`, `data/perturbseq/pu/pu_model_matrix.parquet` (`gwas_score` and `gene_burden_score` columns) | Redistributed evidence tables are committed. The opt-in `make vignette-panelc` target re-queries the live API and asserts release 26.06, aborting under any other release. Consult Open Targets directly for reuse terms. |
| HGNC complete gene set (https://www.genenames.org/download/, accessed 2026-07-16) | Gene symbol reference: `symbol`, `status`, `alias_symbol`, `prev_symbol` | `data/perturbseq/sources/hgnc_symbol_subset.txt` | Committed as a reduced column subset of the upstream complete set, not the verbatim upstream release (see the HGNC section below). Release identified by content: 45,021 approved records, latest `date_modified` 2026-07-14. Consult HGNC for reuse terms. |
| STRING (https://string-db.org) | Protein-protein interaction edges (v12.0 protein links, `protein.info`, `protein.aliases`) backing the discordance-network coherence test | Not committed; downloaded to `data/raw/discordance/` by `make discordance-network` (about 83 MB of protein links plus alias/info maps). Derived edge and coherence tables are committed under `figure_data/` (for example `discordance_network.csv`, `discordance_network_edges.csv`). | Raw archives are not redistributed in this repository; they are fetched on demand. Pinned at v12.0. Consult STRING for reuse terms. |
| Gene Ontology (http://geneontology.org) | UniProt-GOA GAF and a UniProt annotation-depth table, used for GO over-representation and dark-proteome statistics (Figure S2) | Not committed; downloaded to `data/raw/discordance/` by `make discordance-network` (about 11 MB GAF, about 3.5 MB annotation-depth table). Derived enrichment tables are committed under `figure_data/`. | Raw files are not redistributed; fetched on demand. Release unpinned (latest available at fetch time); exact reproduction requires the archived copies in `data/raw/discordance/` rather than a re-fetch. Consult UniProt-GOA/EBI for reuse terms. |
| GWAS Catalog / Open Targets Genetics | Per-gene GWAS Catalog hit count (`gwas_total` = 183) and the Open Targets GWAS-evidence disease rollup (`n_diseases_with_gwas_evidence`) | `data/case_vignette/vignette_meta.json` (`gwas_total`, hand-maintained, carried over from a read-only seed file, not re-derivable offline); `figure_data/vignette_panelC_gwas_meta.json` (opt-in output, Open Targets release 26.06) | `gwas_total` has no producing script in this repository. The Open Targets rollup is re-derived by the opt-in `make vignette-panelc` target. Consult the respective providers for reuse terms. |
| Perturb-seq dataset (genome-scale CRISPRi Perturb-seq, primary human CD4+ T cells, Zhu, Dann et al., 2025; https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq) | Knockdown `zscore` signatures (log2FC / lfcSE), PU-model feature matrix inputs, cytokine and cytokine-receptor significance counts | `data/perturbseq/knn/` (gene-name lists committed; the three signature `.npy` matrices are not, see below), `data/perturbseq/pu/pu_model_matrix.parquet`, `data/perturbseq/pu/pu_labels.csv`, `data/perturbseq/sources/cytokine_significant_counts.csv`, `data/perturbseq/sources/cytokine_receptor_significant_counts.csv` | Released via the CZI Virtual Cells Platform at `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad`. Must be cited in any work using these signatures. Consult the CZI Virtual Cells Platform and Zhu, Dann et al. (2025) for reuse terms. |
| T-helper polarization differential expression (https://github.com/emdann/GWT_perturbseq_analysis_2025) | Observational differential-expression z-scores of polarized T-helper subsets against resting Th0 | `data/perturbseq/pu/pu_model_matrix.parquet` (`zscore_Th1`, `zscore_Th2`, `zscore_Th17`, `zscore_Treg` columns) | Redistributed as feature columns in the committed matrix. Consult the repository and Zhu, Dann et al. (2025) for reuse terms. |
| Inborn errors of immunity (IEI) gene list (IUIS, 2024; https://iuis.org/committees/iei/) | Gene-level IEI status, used as a feature and in enrichment tests | `data/reference/IEI_gene_list.csv` | Committed. Consult IUIS for reuse terms. |
| Druggable genome (Finan et al., 2017, https://doi.org/10.1126/scitranslmed.aag1166) | Gene-level druggability annotation | `data/reference/druggable_genome_gene_list.csv` | Committed. Consult the publication for reuse terms. |
| Immune phecode list | 16 immune phecodes hand-curated from the phecodes of the GPS dataset, used to restrict GPS scores to immune indications | `data/reference/immune_phecodes_curated.csv` | Committed. Hand-curated, with no producing script in this repository. |
| gnomAD genic intolerance metrics (v4.1.0, https://gnomad.broadinstitute.org/data) | Missense constraint (`mis.z_score`) and loss-of-function intolerance (`lof.oe_ci.upper`, LOEUF) | `data/perturbseq/pu/pu_model_matrix.parquet` (`mis.z_score` and `lof.oe_ci.upper` columns) | Redistributed as feature columns in the committed matrix. Consult gnomAD for reuse terms. |
| GPS (Genetic Priority Score) dataset (https://zenodo.org/records/10095684) | Per-gene, per-phecode priority scores (Supplementary Table 16 of the GPS paper) | `data/specificity/gps_per_gene.csv`, `data/specificity/gps_drug_labeled_genes.csv` (committed, derived); source workbook `data/raw/gps/GPS_allgenes.xlsx` (not committed) | Derived files are committed. The source workbook (about 5.6 MB) is fetched on demand by `make gps-inputs`. Consult the GPS paper and its authors for reuse terms. |
| Drug trial date records (ChEMBL, https://www.ebi.ac.uk/chembl/; ClinicalTrials.gov, https://clinicaltrials.gov) | First-trial year per drug and approval year per gene, used to construct temporal holdout labels | `data/prospective/drug_trial_dates.json` (keyed by ChEMBL drug ID, with ClinicalTrials.gov NCT identifiers), `data/prospective/positives_dated.csv` (approval years; 332 rows sourced from the ChEMBL REST API, 8 with no recorded source) | Committed. External source with no producing script in this repository. Consult ChEMBL and ClinicalTrials.gov for reuse terms. |

## Data not committed to this repository

### Perturb-seq signature matrices

Three knockdown-signature matrices are float32 memmaps rebuilt from the perturb-seq `.h5ad`
source and are not committed:

| File | Size | Dimensions |
|---|---|---|
| `data/perturbseq/knn/signatures_Rest.npy` | 281 MB | 6,838 x 10,282 |
| `data/perturbseq/knn/signatures_Stim8hr.npy` | 292 MB | 7,112 x 10,282 |
| `data/perturbseq/knn/signatures_Stim48hr.npy` | 294 MB | 7,156 x 10,282 |

Rows are perturbed genes with a significant on-target knockdown; columns are the 10,282 measured
genes; values are the `zscore` layer (log2FC / lfcSE). Obtain them with:

```
make knn-signatures
```

which streams the `zscore` layer of the upstream `.h5ad` file over HTTP range requests from
`s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad` (about 3 to 4 minutes,
network required). Only `signatures_Stim48hr.npy` is read downstream, by
`build_vignette_panelD_knn.py` for panel d of the STAT4 vignette; all three conditions are read by
`build_knn_nearest_target.py` for Supplementary Table 15. The row and column gene-name lists that
accompany the matrices (`signatures_Rest.genes.txt`, `signatures_Stim8hr.genes.txt`,
`signatures_Stim48hr.genes.txt`, `var_gene_names.txt`) are committed, so row and column identity
can be checked without the matrices themselves. Committing the matrices was judged impractical:
Git LFS's 1 GB/month free bandwidth budget would cover only about three clones of the largest file
before throttling.

### `data/raw/`

Everything under `data/raw/` is gitignored (`/data/raw/*`) and absent on a clean clone:

- `data/raw/discordance/`: STRING v12.0 protein links (about 83 MB) plus STRING's
  `protein.info` and `protein.aliases` maps, the UniProt-GOA GAF (about 11 MB), and a UniProt
  annotation-depth table (about 3.5 MB); about 135 MB total. Obtained with
  `make discordance-network` (about 10 minutes on first run, about 1 minute on later runs, since
  downloads are reused and are not removed by `make clean`). These archives back the eight
  Supplementary Table S2 statistics that `make tables` does not rewrite.
- `data/raw/gps/GPS_allgenes.xlsx`: the GPS paper's source workbook (about 5.6 MB). Obtained with
  `make gps-inputs`, which regenerates `data/specificity/gps_per_gene.csv` and
  `data/specificity/gps_drug_labeled_genes.csv` from it.

## The HGNC file specifically

`data/perturbseq/sources/hgnc_symbol_subset.txt` is a column subset of the upstream HGNC complete
set, not the verbatim upstream release. It retains 4 of the upstream file's 54 columns
(`symbol`, `status`, `alias_symbol`, `prev_symbol`) for all 45,021 rows. The subset was taken
because those are the only columns any consumer in this repository reads. Any reuse that needs a
column outside this set requires fetching the complete upstream file from HGNC directly.

HGNC does not version its complete set, so the release is identified by content: **45,021 approved
records, latest `date_modified` 2026-07-14** (md5 `08386692953902f3880b52ba7d8b2f02` for the
54-column upstream file this was projected from). `load_hgnc` prints both on every load, so the
build log records which reference a run harmonized against.

This matters more than a normal provenance note. HGNC renames symbols between releases, and three
genes in this universe are affected: `IL1RL1`, `IL1RL2` and `C17orf99` are approved here but become
`IL33R`, `IL36R` and `IL40` in later releases. `IL1RL2` is an approved drug target, so harmonizing
against a different release silently re-keys it out of every join on `gene`. `load_hgnc` therefore
refuses to fetch the current release unless called with `allow_download=True`; a missing pinned
file raises rather than falling back to the network.

## Licensing

The code in this repository is licensed under the MIT license; see `LICENSE`, with upstream attributions in `NOTICE.md`. Redistributed
third-party data files remain under their original providers' terms, and this repository does not
relicense them. Users intending to reuse or redistribute any external dataset named above
(Open Targets, HGNC, STRING, Gene Ontology / UniProt-GOA, GWAS Catalog, the perturb-seq dataset
of Zhu, Dann et al. (2025), gnomAD, the GPS dataset, or any other source listed here) should
consult that provider's own terms and, where applicable, cite the corresponding publication or
resource independently of this repository.
