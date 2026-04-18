# Single-Cell RNA-seq Analysis Pipeline

> A complete single-cell RNA sequencing (scRNA-seq) analysis project covering raw data preprocessing, downstream clustering analysis, and data structure exploration using industry-standard bioinformatics tools.

---

## Project Overview

| Property | Details |
|----------|---------|
| **Dataset** | 1k & 3k PBMCs (Peripheral Blood Mononuclear Cells) from a Healthy Donor |
| **Source** | 10x Genomics |
| **Tools Used** | Galaxy, STARsolo, DropletUtils, Scanpy, AnnData, Python |
| **Environment** | Galaxy EU · Google Colab |
| **Language** | Python 3.12 |

---

## Repository Structure

```
scrna-seq-analysis/
│
├── galaxy_preprocessing/
│   ├── outputs/               # Count matrices and QC reports from STARsolo
│   ├── workflow/              # Galaxy workflow screenshots
│   └── README.md              # Galaxy preprocessing documentation
│
├── scanpy_analysis/
│   ├── notebooks/             # Full Jupyter/Colab analysis notebook
│   ├── figures/               # All output plots and visualizations
│   └── README.md              # Scanpy analysis documentation
│
├── anndata_exploration/
│   ├── anndata_demo.ipynb     # AnnData tutorial notebook
│   └── README.md              # AnnData structure documentation
│
└── README.md                  ← You are here
```

---

## Section 1 — 10X Preprocessing (Galaxy)

### Overview

Raw **FASTQ files** from a 10X Chromium single-cell experiment were preprocessed using the **Galaxy platform**. The goal was to convert millions of raw sequencing reads into a structured **gene expression count matrix** — the essential starting point for all downstream analysis.

### Steps Performed

1. **Data Upload** — Uploaded paired-end FASTQ files and a 10X cell barcode whitelist to Galaxy
   - *R1*: Contains 16 bp cell barcode + 12 bp UMI
   - *R2*: Contains the actual cDNA sequence (gene read)

2. **Alignment & Quantification** — Ran **RNA STARsolo** to:
   - Align reads to the human reference genome (*hg19*)
   - Match cell barcodes against the 10X whitelist
   - Count UMIs per gene per valid cell barcode

3. **Quality Control** — Used **MultiQC** to inspect:
   - % reads uniquely mapped to genome
   - % reads with valid barcodes
   - Number of cells detected

4. **Cell Filtering** — Applied **DropletUtils** in two ways:
   - *DefaultDrops* method (Cell Ranger equivalent)
   - *EmptyDrops* method (custom threshold filtering)

5. **Output Export** — Downloaded the 3-file count matrix bundle from Galaxy history

### Key Outputs

| File | Format | Description |
|------|--------|-------------|
| `matrix.mtx` | Market Exchange (MTX) | Sparse UMI count matrix |
| `barcodes.tsv` | Tab-separated | One valid cell barcode per line |
| `features.tsv` | Tab-separated | Gene Ensembl ID + gene symbol |
| MultiQC Report | HTML | Interactive alignment QC summary |

---

## Section 2 — Scanpy Analysis

### Overview

The count matrix was loaded into **Python** and analyzed using **Scanpy** — the standard single-cell analysis library. The complete pipeline from raw counts to annotated cell type clusters was performed on the **PBMC 3k dataset**.

### Workflow Summary

```
Raw counts → QC Filtering → Normalization → Log Transform
→ Highly Variable Genes → Regress Out → Scaling → PCA
→ KNN Graph → Leiden Clustering → UMAP → Marker Genes → Cell Type Annotation
```

### Methods Used

| Step | Function | Purpose |
|------|----------|---------|
| **Cell filtering** | `sc.pp.filter_cells` | Remove empty droplets (< 200 genes) |
| **Gene filtering** | `sc.pp.filter_genes` | Remove rarely expressed genes |
| **QC metrics** | `sc.pp.calculate_qc_metrics` | Compute mitochondrial % per cell |
| **Normalization** | `sc.pp.normalize_total` | Correct for sequencing depth |
| **Log transform** | `sc.pp.log1p` | Stabilize variance |
| **HVG selection** | `sc.pp.highly_variable_genes` | Retain ~2000 informative genes |
| **Regression** | `sc.pp.regress_out` | Remove technical confounders |
| **Scaling** | `sc.pp.scale` | Zero mean, unit variance per gene |
| **PCA** | `sc.tl.pca` | Reduce to 40 principal components |
| **KNN graph** | `sc.pp.neighbors` | Build cell similarity graph |
| **Clustering** | `sc.tl.leiden` | Detect cell communities |
| **Visualization** | `sc.tl.umap` | 2D embedding of cell landscape |
| **Marker genes** | `sc.tl.rank_genes_groups` | Identify cluster-defining genes |

### Key Results

- ***2,638 cells*** and ***1,838 genes*** retained after quality control filtering
- ***8 distinct cell clusters*** identified by the Leiden algorithm
- Cell types annotated using known PBMC marker genes:

| Cluster | **Cell Type** | *Key Marker Genes* |
|---------|-------------|------------------|
| 0 | CD4 T cells | *IL7R, CCR7* |
| 1 | CD14 Monocytes | *LYZ, CD14* |
| 2 | B cells | *CD79A, MS4A1* |
| 3 | CD8 T cells | *CD8A, CD8B* |
| 4 | NK cells | *GNLY, NKG7* |
| 5 | CD14 Monocytes | *S100A8, LGALS3* |
| 6 | Dendritic cells | *FCER1A, CST3* |
| 7 | FCGR3A Monocytes | *FCGR3A, MS4A7* |
| 8 | Platelets | *PPBP* |

### Key Figures Produced

- 📊 QC violin plots — genes per cell, total counts, % mitochondrial reads
- 📊 Highly variable genes dispersion plot
- 📊 PCA variance ratio (elbow plot)
- 📊 UMAP colored by Leiden cluster
- 📊 UMAP colored by marker gene expression
- 📊 Dotplot of marker genes across clusters
- 📊 Final annotated UMAP with cell type labels

---

## Section 3 — AnnData Overview

### Overview

**AnnData** (*"Annotated Data"*) is the core data structure used throughout this entire project. It provides a unified, memory-efficient container for the expression matrix and all associated metadata — making the analysis organized and reproducible.

### AnnData Structure

```
AnnData object  (n_obs × n_vars)
│
├── .X          →  Expression matrix  (cells × genes)
│                  Sparse CSR matrix or dense numpy array
│
├── .obs        →  Cell-level metadata  (Pandas DataFrame)
│                  e.g., cluster labels, QC metrics, cell type names
│
├── .var        →  Gene-level metadata  (Pandas DataFrame)
│                  e.g., highly_variable flags, mean expression
│
├── .obsm       →  Multi-dimensional cell embeddings  (dict)
│                  e.g., X_pca  (n_cells × 40)
│                        X_umap (n_cells × 2)
│
├── .obsp       →  Cell-cell graphs  (dict of sparse matrices)
│                  e.g., KNN connectivities, distances
│
├── .layers     →  Alternative expression matrices  (dict)
│                  e.g., "counts" (raw), "log1p" (normalized)
│
└── .uns        →  Unstructured metadata  (dict)
                   e.g., cluster color palettes, analysis parameters
```

### How AnnData Was Used in This Project

| **Slot** | *Contents After Full Analysis* |
|---------|-------------------------------|
| `.X` | Scaled, normalized expression (HVGs only) |
| `.raw` | Pre-HVG normalized counts (preserved for DE analysis) |
| `.obs` | `n_genes_by_counts`, `pct_counts_mt`, `leiden`, `cell_type` |
| `.var` | `mt` flag, `highly_variable` flag, `means`, `dispersions` |
| `.obsm["X_pca"]` | 40-dimensional PCA embedding |
| `.obsm["X_umap"]` | 2-dimensional UMAP coordinates |
| `.obsp["connectivities"]` | KNN graph (sparse, n_cells × n_cells) |
| `.uns["leiden_colors"]` | Color palette for 8 Leiden clusters |
| `.uns["rank_genes_groups"]` | Wilcoxon marker gene test results |

---

## 🔧 Requirements

```bash
pip install scanpy anndata matplotlib seaborn pandas numpy scipy
pip install python-igraph leidenalg
```

---

## 📖 References

- [Galaxy Training — Pre-processing of 10X Single-Cell RNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
- [Scanpy PBMC3k Tutorial](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html)
- [AnnData Getting Started](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [scverse AnnData Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)
- *Wolf et al. (2018), Scanpy: large-scale single-cell gene expression data analysis, Genome Biology*
- *Tekman et al. (2020), A single-cell RNA-sequencing training and analysis suite using the Galaxy framework, GigaScience*

---

## 👤 Author

**[Your Name]**
*Bioinformatics Assignment — Single-Cell RNA-seq Analysis*
*April 2026*
