# Part 2: Basic scRNA-seq Preprocessing and Clustering

**Notebook:** `basic-scrna-tutorial_updated.ipynb`  
**Original author:** [Taha Shmi](https://github.com/tahashmi/algos/blob/main/basic-scrna-tutorial_updated.ipynb), adapted from the [scverse tutorials](https://scverse-tutorials.readthedocs.io/en/latest/)  
**Based on:** [OpenProblems NeurIPS 2021 benchmarking dataset](https://openproblems.bio/competitions/neurips_2021/) — bone marrow mononuclear cells from healthy human donors, measured with 10X Multiome

---

## Overview

This notebook implements a complete single-cell RNA-seq analysis pipeline starting from a raw count matrix and ending with annotated cell types. It uses [Scanpy](https://scanpy.readthedocs.io/) as the core analysis framework and covers every standard step: quality control, normalisation, dimensionality reduction, clustering, and cell type annotation. The dataset contains **8,785 cells** across two donors (`s1d1`, `s1d3`) with **36,601 genes** measured.

---

## Pipeline

### 1. Data Loading
Two samples are downloaded from the scverse example data repository as `.h5` files and read in with `sc.read_10x_h5`. They are concatenated into a single AnnData object with a `sample` label in `.obs`, and gene names are made unique.

### 2. Quality Control
Three gene categories are flagged in `.var` using name prefixes: mitochondrial (`MT-`), ribosomal (`RPS`/`RPL`), and hemoglobin (`HB`). `sc.pp.calculate_qc_metrics` then computes per-cell statistics including total counts, number of expressed genes, and percentage of counts from each gene group.

**QC outputs:**
- Violin plots of `n_genes_by_counts`, `total_counts`, and `pct_counts_mt`
- Scatter plot of total counts vs. gene count, coloured by mitochondrial percentage

A permissive filter is then applied: cells with fewer than 100 expressed genes and genes detected in fewer than 3 cells are removed. This conservative approach preserves cells that might otherwise be removed prematurely — refinement happens later after clustering.

### 3. Doublet Detection
Scrublet (`sc.external.pp.scrublet`) simulates artificial doublets and scores each real cell by how similar it is to them. The result is a `doublet_score` and a binary `predicted_doublet` flag added to `.obs`. Doublets are not removed immediately — instead they are visualised on the UMAP after clustering so that entire clusters with elevated doublet scores can be identified and removed more robustly (see Section 7).

### 4. Normalisation
Raw counts are saved to `adata.layers["counts"]` before any transformation. `sc.pp.normalize_total` then scales each cell to the same total count depth (median normalisation), and `sc.pp.log1p` applies a log(x+1) transformation. This produces values in `.X` that are comparable across cells and approximately normally distributed — required for PCA and downstream statistics.

### 5. Feature Selection
`sc.pp.highly_variable_genes` identifies the 2,000 most variable genes across the dataset, computed per sample (using `batch_key="sample"`) to avoid selecting genes that are merely variable due to batch effects. Only these genes are used in PCA.

**Output:** A plot showing mean expression vs. dispersion, with the selected highly variable genes highlighted.

### 6. Dimensionality Reduction
PCA is run on the highly variable genes. A variance ratio plot (`sc.pl.pca_variance_ratio`) is used to inspect how many PCs capture meaningful variation — this informs the number of PCs used in the neighbourhood graph step.

A kNN neighbourhood graph is then computed from the PCA embedding (`sc.pp.neighbors`), and UMAP coordinates are derived from this graph (`sc.tl.umap`).

**Output:** UMAP coloured by `sample` — used to assess batch effects. In this dataset only a minor batch effect is observed, so no integration is performed.

### 7. Clustering
Leiden clustering is applied at three resolutions (0.02, 0.5, 2.0), producing three alternative cluster assignments stored in `.obs`. These are visualised side-by-side on the UMAP to choose an appropriate resolution — too low merges distinct populations, too high over-fragments them. Resolution 0.5 is selected as the working clustering.

After clustering, predicted doublets are overlaid on the UMAP. Cells flagged as doublets are then removed with a boolean mask, and QC metrics (total counts, mitochondrial %, gene count) are re-visualised to confirm the remaining cells look healthy.

### 8. Cell Type Annotation

Three complementary approaches are used:

**Manual marker genes** — A curated dictionary of 16 immune cell types and their marker genes is defined (e.g. CD14+ Monocytes: `FCN1`, `CD14`; NK cells: `GNLY`, `NKG7`). `sc.pl.dotplot` visualises the expression of all markers across clusters, allowing manual inspection of which cluster corresponds to which cell type.

**CellTypist (automatic)** — A pre-trained logistic regression model (`Immune_All_Low.pkl`) trained on 20 tissues from 18 studies is applied. With `majority_voting=True`, each Leiden cluster receives the cell type label voted by the majority of its constituent cells. The result (`majority_voting`) is stored in `.obs` and visualised on the UMAP.

**Decoupler + PanglaoDB (enrichment-based)** — Canonical marker genes are retrieved from PanglaoDB via decoupler and filtered to human genes. Multivariate linear regression (`dc.run_mlm`) scores each cell against each cell type's marker programme. The resulting t-value scores are stored in `.obsm["mlm_estimate"]` and visualised as a UMAP coloured by enrichment score for each major cell type. Each cluster is then labelled by the cell type with the highest mean score, producing a `dc_anno` column in `.obs`.

**Output:** A final UMAP comparing `majority_voting` (CellTypist) and `dc_anno` (decoupler) side-by-side — the two annotations are largely concordant.

### 9. Differentially Expressed Marker Genes
`sc.tl.rank_genes_groups` runs Wilcoxon tests for each cluster against all others, identifying genes that are significantly upregulated per cluster. Results are filtered to fold-change ≥ 1.5. The top 5 DE genes per cluster are shown in a dot plot.

A specific cluster (`leiden_res0_5` cluster 3) is highlighted: its top markers (`LYZ`, `ACTB`, `S100A6`, `S100A4`, `CST3`) are visualised on the UMAP and in violin plots across all clusters — demonstrating how DE genes can be used to validate or refine cell type calls.

---

## Repository Structure

```
part2-scrna-scanpy/
├── README.md
└── basic-scrna-tutorial_updated.ipynb   ← This notebook
```

---

## Dependencies

```bash
pip install scanpy anndata pooch scrublet leidenalg celltypist decoupler omnipath
```

