# Scanpy scRNA-seq Analysis

> Complete downstream single-cell RNA-seq analysis of the **PBMC 3k dataset** using **Scanpy** in Python — from raw count loading through to fully annotated cell type clusters.

---

## Goal

***Transform a raw count matrix → biologically meaningful, annotated cell type clusters with visualizations***

---

## Folder Contents

```
scanpy_analysis/
│
├── notebooks/
│   └── pbmc3k_analysis.ipynb       # Full Scanpy analysis pipeline (Google Colab)
│
├── figures/
│   ├── 01_qc_violin.png            # QC metrics — genes, counts, % mitochondrial
│   ├── 02_qc_scatter.png           # Total counts vs % mitochondrial reads
│   ├── 03_highly_variable.png      # Highly variable gene dispersion plot
│   ├── 04_pca_variance.png         # PCA variance ratio (elbow plot)
│   ├── 05_umap_markers.png         # UMAP colored by marker gene expression
│   ├── 06_umap_leiden.png          # UMAP colored by Leiden cluster
│   ├── 07_dotplot.png              # Marker genes dotplot across clusters
│   └── 08_umap_annotated.png       # Final UMAP with cell type labels
│
└── README.md
```

---

## Dataset

| Property | Details |
|----------|---------|
| **Sample** | 3k PBMCs from a Healthy Donor |
| **Full name** | Peripheral Blood Mononuclear Cells |
| **Source** | `sc.datasets.pbmc3k()` — 10x Genomics |
| **Cells (raw)** | 2,700 |
| **Genes (raw)** | 32,738 |
| **Cells after QC** | ***2,638*** |
| **Genes after QC** | ***~1,838*** |

---

## Analysis Pipeline

```
sc.datasets.pbmc3k()          ← Load dataset
        │
        ▼
Filter cells & genes          ← Remove low quality cells and rare genes
        │
        ▼
QC metrics & thresholds       ← Flag and remove dying cells (high % mito)
        │
        ▼
Normalize → Log1p transform   ← Make cells comparable
        │
        ▼
Highly variable genes         ← Keep only the most informative genes
        │
        ▼
Regress out confounders       ← Remove technical variation
        │
        ▼
Scale                         ← Zero mean, unit variance
        │
        ▼
PCA (40 components)           ← Compress into key dimensions
        │
        ▼
KNN neighbor graph            ← Connect similar cells
        │
        ▼
Leiden clustering             ← Find cell communities
        │
        ▼
UMAP visualization            ← 2D view of cell landscape
        │
        ▼
Rank genes (Wilcoxon)         ← Find marker genes per cluster
        │
        ▼
Cell type annotation          ← Label clusters biologically
```

---

## Methods

### **Quality Control**

Three filters were applied to remove low quality cells:

| Filter | Threshold | Removes |
|--------|-----------|---------|
| Minimum genes per cell | > 200 | *Empty droplets* |
| Maximum genes per cell | < 2,500 | *Likely doublets* |
| Mitochondrial read % | < 5% | *Dying/stressed cells* |

> ***Why mitochondrial reads?*** When a cell is dying, cytoplasmic RNA leaks out but mitochondrial RNA (enclosed in mitochondria) stays behind. A high % of mitochondrial reads therefore indicates a low quality or dead cell.

---

### **Normalization**

Library-size normalization was applied by scaling each cell to **10,000 total counts**, followed by ***log(x+1) transformation*** to compress the dynamic range and stabilize variance across genes.

---

### **Highly Variable Gene (HVG) Selection**

Genes were ranked by their **mean expression** and **dispersion** (variance relative to mean). Only genes with:
- `min_mean = 0.0125`, `max_mean = 3`
- `min_disp = 0.5`

were retained, reducing the feature space from ~32,000 genes to ***~2,000 highly informative genes***.

---

### **Regression & Scaling**

Technical confounders (`total_counts` and `pct_counts_mt`) were regressed out using `sc.pp.regress_out()`. Each gene was then scaled to **zero mean and unit variance**, clipped at a maximum value of 10 to prevent outlier dominance.

---

### **Dimensionality Reduction**

**PCA** was applied to the scaled HVG matrix using 40 principal components (`svd_solver='arpack'`). The PCA variance ratio plot was used to confirm that 40 PCs captured the majority of variance.

**UMAP** was computed from the KNN neighbor graph (k=10, n_pcs=40), projecting cells into 2D while preserving local structure.

---

### **Clustering**

The **Leiden algorithm** was applied at `resolution=0.9`, yielding ***8 distinct clusters***. Higher resolution values produce more clusters; lower values produce fewer. The Leiden algorithm improves upon the older Louvain algorithm by guaranteeing well-connected communities.

---

### **Marker Gene Identification**

Differentially expressed genes were identified for each cluster using the **Wilcoxon rank-sum test** via `sc.tl.rank_genes_groups`. The top 25 marker genes per cluster were visualized.

---

## Results

### Cell Type Annotations

| Cluster | **Cell Type** | *Key Marker Genes* | Biology |
|---------|-------------|------------------|---------|
| 0 | CD4 T cells | *IL7R, CCR7* | Helper T cells |
| 1 | CD14 Monocytes | *LYZ, CD14* | Innate immune cells |
| 2 | B cells | *CD79A, MS4A1* | Antibody-producing cells |
| 3 | CD8 T cells | *CD8A, CD8B* | Cytotoxic T cells |
| 4 | NK cells | *GNLY, NKG7* | Natural killer cells |
| 5 | CD14 Monocytes | *S100A8, LGALS3* | Inflammatory monocytes |
| 6 | Dendritic cells | *FCER1A, CST3* | Antigen-presenting cells |
| 7 | FCGR3A Monocytes | *FCGR3A, MS4A7* | Patrolling monocytes |
| 8 | Platelets | *PPBP* | Blood clotting cells |

---

## How to Reproduce

### Install dependencies

```bash
pip install scanpy anndata matplotlib seaborn
pip install python-igraph leidenalg
```

### Run in Google Colab

1. Open `notebooks/pbmc3k_analysis.ipynb` in Google Colab
2. Run all cells in order (Runtime → Run all)
3. Figures will be saved to the `figures/` folder

### Key code snippet

```python
import scanpy as sc

# Load dataset
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

# QC filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize and cluster
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ... see full notebook for complete pipeline
```

---

## 🔗 References

- [Scanpy Documentation](https://scanpy.readthedocs.io)
- [Scanpy PBMC3k Tutorial](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html)
- *Wolf et al. (2018), Scanpy: large-scale single-cell gene expression data analysis, Genome Biology*
- *Satija et al. (2015), Spatial reconstruction of single-cell gene expression data, Nature Biotechnology*
