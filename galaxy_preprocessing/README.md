# Galaxy 10X Preprocessing

> Preprocessing of raw 10X Genomics scRNA-seq FASTQ files into a structured gene expression count matrix using the **Galaxy platform** and **STARsolo**.

---

## Goal

***Convert raw sequencing reads → clean, filtered count matrix (cells × genes)***

---

## Folder Contents

```
galaxy_preprocessing/
│
├── outputs/
│   ├── matrix.mtx              # Sparse UMI count matrix (STARsolo output)
│   ├── barcodes.tsv            # Valid cell barcodes detected
│   ├── features.tsv            # Gene list (Ensembl ID + symbol)
│   └── multiqc_report.html     # Interactive QC summary report
│
├── workflow/
│   └── starsolo_workflow.png   # Screenshot of Galaxy workflow
│
└── README.md
```

---

## 🧬 Dataset

| Property | Details |
|----------|---------|
| **Sample** | 1k PBMCs from a Healthy Donor |
| **Full name** | Peripheral Blood Mononuclear Cells |
| **Source** | 10x Genomics (Zenodo record 3457880) |
| **Chemistry** | 10X Chromium ***v3*** |
| **Reference genome** | *hg19 (GRCh37)* |
| **Input files** | 4 FASTQ files across 2 sequencing lanes (L001, L002) |

---

## 🧪 10X Chemistry — Understanding the Input Files

The 10X Chromium v3 system produces **3 reads per lane**:

| File | Read | Contains | Length |
|------|------|----------|--------|
| `R1` | Read 1 | **Cell Barcode** (16 bp) + **UMI** (12 bp) | 28 bp |
| `R2` | Read 2 | cDNA sequence (actual gene read) | 91 bp |
| `I1` | Index | Sample index (for demultiplexing lanes) | 8 bp |

> ***Note:*** STARsolo only requires R1 and R2 — the I1 index file is not needed.

### Key Concepts

- **Cell Barcode** — A unique 16-base sequence that identifies *which cell* a read came from
- **UMI (Unique Molecular Identifier)** — A 12-base sequence tagging each RNA molecule *before* PCR amplification, used to remove duplicate reads
- **Whitelist** — A predefined list of ~3.7 million known valid 10X barcodes. Only reads matching this list are counted as real cells

---

## Workflow Steps

### Step 1 — Data Upload & Organization

- Created a new Galaxy history: *"scRNA-seq 10X dataset tutorial"*
- Imported the **sub-sampled FASTQ files** from Zenodo:
  - `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz`
  - `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz`
  - `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz`
  - `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz`
- Imported supporting files:
  - `Homo_sapiens.GRCh37.75.gtf` — Gene annotation file
  - `3M-february-2018.txt.gz` — 10X cell barcode whitelist

---

### Step 2 — Alignment & Quantification with STARsolo

**Tool:** *RNA STARsolo (Galaxy version 2.7.11a+galaxy1)*

STARsolo simultaneously:
- Aligns **R2** reads to the human reference genome
- Matches **R1** barcodes against the whitelist
- Counts **UMIs** per gene per valid cell barcode

**Key parameters used:**

| Parameter | Value |
|-----------|-------|
| Reference genome | Human (*hg19*) |
| Gene annotation | `Homo_sapiens.GRCh37.75.gtf` |
| Chemistry | Chromium v3 |
| Barcode whitelist | `3M-february-2018.txt.gz` |
| UMI deduplication | CellRanger2-4 algorithm |
| Cell filtering | *Disabled* (filtered manually in next step) |

**STARsolo produced 6 output files:**
- Log file
- Feature Statistic Summaries
- BAM alignment file
- `matrix.mtx` — count matrix
- `barcodes.tsv` — detected cell barcodes
- `genes.tsv` — gene list

---

### Step 3 — Quality Control with MultiQC

**Tool:** *MultiQC (Galaxy version 1.27+galaxy0)*

Assessed alignment quality by inspecting the STARsolo log:

| Metric | Meaning |
|--------|---------|
| `yesWLmatchExact` | Reads with exactly matching barcodes |
| `yesCellBarcodes` | Total cells detected (~5,200 before filtering) |
| `noNoFeature` | Reads that mapped to genome but not to any gene |
| `yessubWLmatch_UniqueFeature` | Reads counted towards a unique gene |

---

### Step 4 — Cell Filtering with DropletUtils

**Tool:** *DropletUtils (Galaxy version 1.10.0+galaxy2)*

Two filtering methods were applied:

#### ***Method A — DefaultDrops (Cell Ranger equivalent)***
- Expected cells: 3,000
- Upper quantile: 0.99
- Lower proportion: 0.1
- Result: ~***272 high quality cells*** recovered

#### ***Method B — EmptyDrops (Introspective/Custom)***
- Barcode rank plot used to find *knee* and *inflection* thresholds
- Lower-bound threshold: 200 UMIs
- FDR threshold: 0.01 (max 1% false positives)
- Result: ~***279 high quality cells*** recovered

> ***Key insight:*** The EmptyDrops method recovers slightly more cells by using a statistical model rather than a fixed threshold, which can improve resolution in downstream clustering.

---

### Step 5 — Output Export

Downloaded the final filtered 3-file bundle from Galaxy history for use in Scanpy:
- `matrix.mtx`
- `barcodes.tsv`
- `features.tsv`

---

## Key Outputs

| File | Format | Description |
|------|--------|-------------|
| `matrix.mtx` | Market Exchange (MTX) | Sparse count matrix — rows=genes, cols=cells |
| `barcodes.tsv` | Tab-separated | One valid cell barcode per line |
| `features.tsv` | Tab-separated | Ensembl gene ID + gene symbol |
| `multiqc_report.html` | HTML | Interactive alignment QC report |

---

## Why MTX Format?

The count matrix is stored as a **sparse matrix** because:
- Most genes are *not expressed* in any given cell
- A full 60,000 genes × 3,000,000 barcodes matrix would be enormous
- Sparse format stores ***only the non-zero values***, saving significant memory

---

## 🔗 References

- [Galaxy Tutorial — Pre-processing of 10X Single-Cell RNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
- *Tekman et al. (2020), A single-cell RNA-sequencing training and analysis suite using the Galaxy framework, GigaScience*
- [Zenodo Dataset Record 3457880](https://zenodo.org/record/3457880)

