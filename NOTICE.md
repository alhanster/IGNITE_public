# Notices

The MIT license in `LICENSE` covers this repository's code. It does not extend to the upstream
data described below, which remains under each source's own terms. `DATA_AVAILABILITY.md` is the
authoritative per-source record of what is committed, what is fetched on demand, and the terms
each is used under.

## Upstream scientific data source

The perturb-seq feature layers consumed by this repository (the committed
`data/perturbseq/pu/pu_model_matrix.parquet`, and the knockdown signature matrices rebuilt by
`src/analysis/vignette/build_knn_signatures.py`) derive from the genome-scale CRISPRi Perturb-seq resource
in primary human CD4+ T cells of Zhu, Dann et al. (2025), released via the CZI Virtual Cells
Platform:

    s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad

That resource is the upstream scientific data source for every functional-genomic feature in this
model and must be cited in any work using them. It remains under its own terms.

## Attribution

Portions of the perturb-seq feature-engineering approach are derived from
emdann/GWT_perturbseq_analysis_2025 (MIT License, (c) 2025 Emma Dann). The scripts that implement
it are not part of this reproduction subset (only their committed output,
`pu_model_matrix.parquet`, is), but the attribution travels with the data.

## Redistributed third-party data

Files under `data/` include reference tables redistributed from Open Targets, HGNC, the GWAS
Catalog, the Gene Ontology and gnomAD, and tables under `figure_data/` include values derived
from STRING. Raw STRING and UniProt-GOA archives are not redistributed; they are fetched on demand.
Each source remains under its own terms and should be cited independently. See
`DATA_AVAILABILITY.md` and `figure_data/PROVENANCE.md`.
