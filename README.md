# scRNA-seq Analysis Pipeline

## Introduction

This repository presents a structured workflow for the analysis of single-cell RNA sequencing (scRNA-seq) data, covering three key stages: preprocessing, data structure understanding, and downstream analysis. For a detailed breakdown of each part, check the respective directories for documentation.

The project is divided into three main parts:

1. **Preprocessing of 10X scRNA-seq data (Galaxy):**
   Raw sequencing reads from a 10X Genomics dataset are processed using the Galaxy platform. This includes alignment and quantification using STARsolo, followed by cell filtering with DropletUtils to generate a high-quality gene-cell count matrix.

2. **AnnData Tutorial (Data Structure Exploration):**
   This section focuses on understanding the AnnData object, which is the foundational data structure used in single-cell analysis workflows. It demonstrates how gene expression data and metadata are stored, accessed, and manipulated.

3. **Scanpy Analysis (Downstream Analysis):**
   Using Scanpy, the processed dataset is analyzed through standard single-cell workflows including quality control, normalization, dimensionality reduction, clustering, and visualization.

---

## Directory Structure and Navigation

The repository is organized into three main directories corresponding to each part of the workflow:

```
scRNAseq_analysis/
├── 01_preprocessing_galaxy
│   ├── Documentation_Galaxy.md
│   ├── Workflow.png
│   ├── data
│   ├── logs
│   └── qc_reports
│
├── 02_anndata_tutorial
│   └── anndata_exploration.ipynb
│
├── 03_scanpy_analysis
│   └── basic-scrna-tutorial_updated.ipynb
│
└── README.md
```

* **01_preprocessing_galaxy/**
  Contains outputs from the Galaxy-based preprocessing pipeline, including the final filtered count matrix (`data/`), quality control reports (`qc_reports/`), and STARsolo logs (`logs/`). The workflow image and documentation file describe the steps performed.

* **02_anndata_tutorial/**
  Includes a Jupyter notebook demonstrating the structure and manipulation of AnnData objects.

* **03_scanpy_analysis/**
  Contains the Scanpy notebook used for downstream analysis and visualization of single-cell data.

---

## Learning Outcomes

This project provides a comprehensive understanding of the scRNA-seq analysis pipeline:

* **Data Preprocessing:**
  Learned how raw sequencing reads are transformed into a structured gene-cell count matrix using alignment and filtering tools. The importance of removing empty droplets and low-quality cells to ensure reliable downstream analysis was emphasized.

* **Understanding Data Structures:**
  Gained familiarity with the AnnData format, including how expression matrices and metadata are stored and accessed. This is critical for efficient manipulation and analysis of large-scale single-cell datasets.

* **Downstream Analysis:**
  Developed practical experience with Scanpy for performing normalization, dimensionality reduction (PCA, UMAP), clustering, and visualization. This highlights how biological insights are extracted from processed data.

* **Reproducibility and Organization:**
  Learned how to structure a bioinformatics project in a clear and reproducible manner using GitHub, including proper documentation and separation of workflow stages.

---

## References

* Galaxy Training Network. *Pre-processing of 10X Single-Cell RNA Datasets*:
  https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html

* Scanpy Tutorials (scverse):
  https://github.com/scverse/scanpy-tutorials

* Updated Scanpy Tutorial:
  https://github.com/tahashmi/algos/blob/main/basic-scrna-tutorial_updated.ipynb

* AnnData Documentation:
  https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

* scverse AnnData Tutorial:
  https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html

* 10X Genomics Dataset (PBMC 1k, v3 chemistry):
  https://zenodo.org/record/3457880

