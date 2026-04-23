# Part 3: Getting Started with AnnData

**Notebook:** `anndata_exploration.ipynb`  
**References:**
- [anndata.readthedocs.io — Getting Started](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [scverse-tutorials — AnnData Getting Started](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)

---

## Overview

`AnnData` is the standard data container for single-cell analysis in Python. It holds the expression matrix alongside all associated metadata in one object, keeping everything indexed and aligned. This notebook works through the core components of AnnData using a simulated 100-cell × 2000-gene dataset, covering construction, annotation, I/O, and memory-efficient access patterns.

---

## Notebook Sections

### 1. Initializing an AnnData Object
Constructs an AnnData object from a sparse count matrix (`scipy.sparse.csr_matrix`). Sparsity is the natural choice for scRNA-seq data since most genes are unexpressed in most cells. The section also covers naming obs/var axes and demonstrates subsetting by name, integer index, and the view vs. copy distinction — important because subsetting returns a read-only view by default, and `.copy()` is needed before any modification.

### 2. Aligned Metadata — `obs` and `var`
`.obs` and `.var` are Pandas DataFrames that stay index-aligned to the cell and gene axes respectively. The notebook adds cell type, batch, and total UMI count to `obs`, and mitochondrial flag, mean expression, and chromosome to `var`. Boolean masking directly on these DataFrames is shown as the standard way to filter cells or genes — e.g. keeping only T cells or dropping mitochondrial genes.

### 3. Observation/Variable-Level Matrices — `obsm`, `varm`, `obsp`
These slots hold multi-dimensional arrays aligned to one axis. `obsm` stores per-cell arrays (PCA coordinates, UMAP embeddings), `varm` stores per-gene arrays (PCA loadings), and `obsp` stores pairwise cell-cell matrices (kNN connectivity graph). Understanding these is necessary because Scanpy writes its results — PCA, UMAP, neighbours — directly into these slots rather than returning separate objects.

### 4. Unstructured Metadata — `uns`
`.uns` is a free-form dictionary for anything that doesn't fit an aligned structure: colour palettes, tool parameters, dataset provenance, clustering run configs. The notebook stores colour maps, neighbour graph parameters, and PCA variance ratios — all examples of what Scanpy populates automatically during a real analysis run.

### 5. Layers
`.layers` stores alternative expression matrices of the same shape as `.X`. The notebook implements the standard preprocessing sequence — saving raw counts, then library-size normalisation (scaling to 10,000 counts per cell), then log1p transformation — as separate named layers. This pattern matters because Scanpy overwrites `.X` at each step, and layers are the correct way to retain earlier representations without duplicating the whole object.

### 6. Conversion to DataFrames
`adata.to_df()` converts `.X` or a named layer to a dense Pandas DataFrame for use with standard data science tooling. The section shows extracting expression values and joining them with `obs` metadata using `groupby` — useful for quick summaries but not scalable to very large datasets due to densification.

### 7. Writing and Reading `.h5ad` Files
`.h5ad` is the standard binary format for AnnData, storing all components in a single HDF5 file. The notebook writes and re-reads a complete object, verifies that all components survive the round-trip (including `.obsm`, `.layers`, `.uns`), and compares uncompressed vs. gzip-compressed file sizes. This format is important because it is the interchange standard between Python and R (via `zellkonverter`) in single-cell workflows.

### 8. Partial Reading of Large Data (Backed Mode)
For datasets that exceed available RAM, `read_h5ad(..., backed="r")` keeps `.X` on disk and only loads metadata into memory. The notebook demonstrates filtering cells by cell type from metadata alone (without loading `.X`), then using `.to_memory()` to pull only the matching subset into RAM. This is the correct pattern for working with production-scale datasets of millions of cells.

### 9. Views vs. Copies
A brief but practically important section: subsetting returns a view — lightweight and read-only. Attempting to modify a view raises a warning. The explicit rule is to always call `.copy()` before adding or changing anything on a subset.

