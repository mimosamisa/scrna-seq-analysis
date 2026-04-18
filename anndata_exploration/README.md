# AnnData Structure & Exploration

> A practical introduction to **AnnData** (*Annotated Data*) — the core data structure powering the entire scverse ecosystem, used throughout this project to store, access, and manipulate single-cell data.

---

## Goal

***Understand how AnnData organizes single-cell data and demonstrate each of its key components using real and simulated data***

---

## Folder Contents

```
anndata_exploration/
│
├── anndata_demo.ipynb      # Standalone AnnData tutorial notebook (Google Colab)
└── README.md
```

---

## What is AnnData?

**AnnData** is a Python data structure specifically designed for ***matrix-like biological data***. It stores the expression matrix alongside all cell and gene metadata in a ***single, well-organized object*** — enabling reproducible, efficient single-cell analysis.

It is the backbone of the entire **scverse ecosystem** — including *Scanpy*, *scVI*, *CellRank*, and more.

> ***Think of it as:*** A smarter version of a spreadsheet, where the rows (cells) and columns (genes) each carry their own metadata, and the whole object can be saved and shared as a single compressed file.

---

## AnnData Structure

```
AnnData object  (n_obs × n_vars)
│
├── .X          →  Core expression matrix  (cells × genes)
│                  Stored as a sparse CSR matrix (memory efficient)
│
├── .obs        →  Cell-level metadata  (Pandas DataFrame)
│                  One row per cell
│                  e.g., cluster labels, QC metrics, cell type names
│
├── .var        →  Gene-level metadata  (Pandas DataFrame)
│                  One row per gene
│                  e.g., highly_variable flag, mean expression
│
├── .obsm       →  Multi-dimensional cell embeddings  (dict of arrays)
│                  e.g., X_pca  →  shape: (n_cells × 40)
│                        X_umap →  shape: (n_cells × 2)
│
├── .varm       →  Multi-dimensional gene matrices  (dict of arrays)
│                  e.g., PCs → PCA gene loadings
│
├── .obsp       →  Cell-cell graphs  (dict of sparse matrices)
│                  e.g., connectivities, distances  (KNN graph)
│
├── .layers     →  Alternative expression matrices  (dict)
│                  Same shape as .X — stores different data versions
│                  e.g., "raw_counts", "log1p"
│
└── .uns        →  Unstructured metadata  (plain Python dict)
                   e.g., cluster color palettes, analysis parameters
```

---

## Slot-by-Slot Reference

### `.X` — *The Expression Matrix*

The **primary data matrix**. Rows are cells (observations), columns are genes (variables). Stored as a ***sparse CSR matrix*** to save memory, since the vast majority of gene counts in any cell are zero.

```python
print(adata.X)                       # sparse matrix summary
print(adata.X[:5, :5].toarray())     # view as dense array
```

---

### `.obs` — *Cell-Level Metadata*

A **Pandas DataFrame** with one row per cell. Holds any per-cell annotation added during analysis. Categorical columns are preferred for efficiency.

```python
adata.obs["cell_type"] = pd.Categorical(["B", "T", "Monocyte", ...])
adata.obs["leiden"] = ...            # cluster labels from Leiden algorithm
print(adata.obs.head())
```

---

### `.var` — *Gene-Level Metadata*

A **Pandas DataFrame** with one row per gene. Holds any per-gene annotation.

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")   # mitochondrial flag
adata.var["highly_variable"]                               # set by HVG step
print(adata.var.head())
```

---

### `.obsm` — *Cell Embeddings*

A **dictionary of 2D arrays**, one row per cell. Used to store dimensionality reduction results.

```python
adata.obsm["X_pca"]    # shape: (n_cells, 40)  → PCA coordinates
adata.obsm["X_umap"]   # shape: (n_cells, 2)   → UMAP coordinates
```

---

### `.layers` — *Alternative Matrices*

A **dictionary of matrices** with the same shape as `.X`. Useful for storing multiple versions of the data side by side without losing any information.

```python
adata.layers["raw_counts"] = adata.X.copy()    # before normalization
adata.layers["log1p"] = normalized_matrix       # after log transform
```

---

### `.uns` — *Unstructured Metadata*

A **plain Python dictionary**. Stores any metadata that doesn't fit the structured slots above — completely flexible.

```python
adata.uns["experiment"] = "PBMC Tutorial"
adata.uns["leiden_colors"]     # cluster color palette (auto-set by Scanpy)
adata.uns["rank_genes_groups"] # Wilcoxon marker gene results
```

---

## How AnnData Was Used in This Project

After completing the full Scanpy pipeline, the PBMC AnnData object contained:

| **Slot** | *Contents* | Added By |
|---------|-----------|----------|
| `.X` | Scaled, normalized expression (HVGs only) | `sc.pp.scale` |
| `.raw` | Pre-HVG normalized counts | `adata.raw = adata` |
| `.obs["n_genes_by_counts"]` | Genes detected per cell | `sc.pp.calculate_qc_metrics` |
| `.obs["pct_counts_mt"]` | % mitochondrial reads per cell | `sc.pp.calculate_qc_metrics` |
| `.obs["leiden"]` | Cluster label (0–8) | `sc.tl.leiden` |
| `.obs["cell_type"]` | Annotated cell type name | *Manual annotation* |
| `.var["mt"]` | Mitochondrial gene flag | *Manual* |
| `.var["highly_variable"]` | HVG selection flag | `sc.pp.highly_variable_genes` |
| `.obsm["X_pca"]` | 40-dim PCA coordinates | `sc.tl.pca` |
| `.obsm["X_umap"]` | 2-dim UMAP coordinates | `sc.tl.umap` |
| `.obsp["connectivities"]` | KNN graph (sparse) | `sc.pp.neighbors` |
| `.uns["leiden_colors"]` | Color palette for 8 clusters | `sc.tl.leiden` |
| `.uns["rank_genes_groups"]` | Wilcoxon test results | `sc.tl.rank_genes_groups` |

---

## Subsetting AnnData

Subsetting follows the same rules as Pandas DataFrames — using names, indices, or boolean masks.

```python
# Subset by cell name
adata[["Cell_1", "Cell_2"], :]

# Subset by boolean mask — only B cells
b_cells = adata[adata.obs["cell_type"] == "B cells", :]

# Subset by cluster
cluster0 = adata[adata.obs["leiden"] == "0", :]
```

> ***Important:*** Subsetting returns a ***view*** (not a copy). Use `.copy()` to create an independent object that can be modified safely.

---

## Saving & Loading `.h5ad`

AnnData uses the **`.h5ad`** file format — a compressed *HDF5* file that stores ***all slots*** (matrix, metadata, embeddings, graphs) together in one file.

```python
# Save everything to disk
adata.write_h5ad("pbmc3k_processed.h5ad")

# Load it back in any session
import anndata as ad
adata = ad.io.read_h5ad("pbmc3k_processed.h5ad")
```

---

## How to Reproduce

```bash
pip install anndata numpy pandas scipy
```

Open `anndata_demo.ipynb` in Google Colab and run all cells in order.

---

## 🔗 References

- [AnnData Documentation](https://anndata.readthedocs.io/en/latest/)
- [Getting Started with AnnData](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [scverse AnnData Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)
- *Virshup et al. (2021), The scverse project provides a computational ecosystem for single-cell omics data analysis, Nature Biotechnology*
