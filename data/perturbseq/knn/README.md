# CD4+ T-cell Perturb-seq knockdown signatures

| File | Size | Committed? | Contents |
|---|---|---|---|
| `signatures_Rest.npy` | 281 MB | **no** | 6,838 × 10,282 `float32` |
| `signatures_Stim8hr.npy` | 292 MB | **no** | 7,112 × 10,282 `float32` |
| `signatures_Stim48hr.npy` | 294 MB | **no** | 7,156 × 10,282 `float32`. Rows = perturbed genes with a significant on-target knockdown; columns = measured genes. Values are the `zscore` layer (log2FC / lfcSE) |
| `signatures_Rest.genes.txt` | 43 KB | yes | Row gene symbols, in row order |
| `signatures_Stim8hr.genes.txt` | 45 KB | yes | Row gene symbols, in row order |
| `signatures_Stim48hr.genes.txt` | 45 KB | yes | Row gene symbols, in row order |
| `var_gene_names.txt` | 65 KB | yes | Column gene symbols, in column order, shared by all three |

All three conditions are written by `make knn-signatures` and all three are read, as described below.

## Not committed, regenerated instead

```bash
make knn-signatures     # OPT-IN, network, ~3-4 min
```

The matrices are 281-294 MB each, two orders of magnitude larger than other files here. Git LFS would serve them from a 1 GB/month free budget, enough for about three clones before it throttles for everyone. `build_knn_signatures.py` streams the `zscore` layer of the source `.h5ad` over HTTP range requests, keeping de-duplicated rows with a significant on-target knockdown (`ontarget_significant`, `n_cells_target >= 50`) per condition.

Nothing downstream depends on the matrices at run time. Both steps that read them, `build_vignette_panelD_knn.py` and `build_knn_nearest_target.py`, have already written their outputs into the committed `figure_data/` contract and skip with a message when the matrices are absent. A clean clone runs `make tables` end to end and passes `make verify-tables` without calling this target. On such a clone, `verify-tables` passes for those five tables without rewriting them.

The four gene-name lists are committed, so row and column identity can be checked without the matrices. Since `make knn-signatures` overwrites them, a clean `git status, data/perturbseq/knn` after a rebuild confirms the source has not moved. This check covers all three `.genes.txt` files, not only the one the vignette needs.

## What reads them, and how far the result is carried

Two steps ask two different questions of the same arrays. Both zero each perturbed gene's own column before L2-normalizing the rows, so a pair's similarity is never inflated by each knockdown being the largest effect in its own row.

- **Panel d of the STAT4 vignette**: `build_vignette_panelD_knn.py`, `Stim48hr` only, since the panel is a stimulated-48-hour analysis. STAT4 against every other perturbed gene, reporting the top six approved immune-drug targets and the 95th percentile of the unlabeled distribution. Produces `figure_data/vignette_panelD_knn.csv` and `vignette_knn_stats.json`.
- **Supplementary Table 15**: `build_knn_nearest_target.py`, then `add_knn_approved_drugs.py`, all three conditions. Every perturbed gene against its nearest approved-target anchor, a maximum over the 105/112/112 anchors that carry a signature in that condition. Produces `figure_data/knn_nearest_target_{Rest,Stim8hr,Stim48hr}.csv`.

The summary statistics from the two tables are not interchangeable. The per-gene numbers agree where they overlap: STAT4's CD3E at 0.355 appears in both, computed independently. Panel d's `knn_p95` is 0.234, while the 95th percentile of Table 15's unlabeled column is 0.354. Both are computed over the same 7,044 unlabeled rows of the same matrix in the same condition, but one is a fixed gene against all genes and the other is each gene against its own best anchor, so the two percentiles measure different quantities.

The cosine is a post-hoc interpretive annotation rather than a screening score: held-out trial-vs-unlabeled discrimination from it is near chance (AUC 0.47-0.53) under every pooling scheme tested. It appears as one panel of a case vignette and one supplementary table, not as a model feature.

## Provenance

Genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells, Zhu, Dann et al. (2025), via the CZI Virtual Cells Platform:

```
s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad
```

33,983 (perturbation x condition) rows by 10,282 measured genes. The resource requires citation in any work using these signatures; see `NOTICE.md`.
