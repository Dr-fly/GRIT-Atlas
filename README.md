# GRIT-Atlas: A Comprehensive Treatment-Induced Resistance Atlas of Glioblastoma

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18009414.svg)](https://zenodo.org/records/18009414)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational pipeline for the **Glioblastoma Resistance Insights from Treatment Atlas (GRIT-Atlas)**. Our study utilizes an integrated single-cell and spatially resolved multi-omic approach to deconstruct the cellular and spatial logic of therapy-induced evolution in IDH-wildtype glioblastoma.

## 🌟 Scientific Highlights
* **Scale and Scope**: We constructed a high-resolution integrated atlas encompassing **978,065 cells** from **296 samples**, capturing the full spectrum of disease evolution from primary diagnosis to recurrence post-standard-of-care (SOC) and combinatorial immunotherapy (ICB + anti-angiogenic therapy).
* **The Spatial Resistance Triad**: Using spatial transcriptomics across 48 patient sections, we identified a core functional unit composed of **cNMF7 (MES-like)** malignant cells, differentiation-arrested **E-MDSCs**, and Type VI Collagen-secreting **myCAFs**.
* **Niche-Specific Resistance**: This triad specifically colonizes the hypoxic microvascular proliferation (**MVP**) and pseudopalisading necrosis (**PAN**) niches.
* **Mechanistic Insight**: We demonstrate that myCAFs function as stromal architects, constructing a fibrotic scaffold via the **Collagen/Fibronectin-CD44 signaling axis** to physically exclude cytotoxic T cells and sustain malignant plasticity.

---

## 📂 Repository Organization

The analysis is modularized into four key computational stages:

| Script Name | Analysis Module | Primary Methodology | Linked Data Files |
| :--- | :--- | :--- | :--- |
| `1_Single_Cell_Preprocessing.ipynb` | Atlas Construction | QC, Covariate Regression, BBKNN | `IntegratedData.zip` |
| `2_Malignant_cNMF_Analysis.ipynb` | Malignant Phenotyping | Consensus NMF (cNMF) | `Malignant.qs` |
| `3_Spatial_Landscape_Analysis.R` | Spatial Modeling & Gradients | RCTD, MistyR, and SPATA2 | `GRIT-Atlas_RCTDres.zip` |
| `4_Niche_Validation.R` | Niche Communications | CellChat Modeling | `CellChat.zip` |

---

## 📊 Data Access (Zenodo)

All processed data objects required to reproduce the findings are hosted on **Zenodo**:
🔗 **[Direct Link: https://zenodo.org/records/18009414](https://zenodo.org/records/18009414)**

### 1. Global Integrated Data
* **`IntegratedData.zip`**: The master Seurat/AnnData object for the complete GRIT-Atlas (978,065 cells). Includes batch-corrected UMAP coordinates, treatment status, clinical response labels, and lineage annotations.

### 2. High-Resolution Lineage Objects
* **`Malignant.qs`**: Details 11 functional states (cNMF), including the resistance-associated **cNMF7 (MES-like)** state and malignant Neuron-like (NEU) states.
* **`Myeloid.zip`**: High-resolution subclustering of 349,583 myeloid cells (23 subpopulations), highlighting **E-MDSC.ADM.HIF1A**.
* **`Lymphocyte.zip`**: Detailed annotations for 36 lymphoid subsets, including cytotoxic T cells and exhausted **Tex.03.ALDOA.MIF**.
* **`Stroma.qs`**: Processed stromal compartment (CAFs, Pericytes, SMCs), specifically annotating the **myCAF.02.COL6** subset.
* **`Endothelial.qs`**: High-resolution mapping of ECs, characterizing immune-resistant subsets (**Vein.VCAM1** and **Arterial.01.DKK2**).

### 3. Spatial Transcriptomics & Communication
* **`GRIT-Atlas_RCTDres.zip`**: RCTD weight matrices and MistyR importance values for 48 glioblastoma patient sections (10x Visium).
* **`CellChat.zip`**: Processed objects documenting communication within the "Spatial Resistance Triad," focused on the **Collagen-CD44 axis**.

---

## 💻 Software Requirements

### **Python Environment (v3.10+)**
* `scanpy` (v1.10.1): QC and normalization.
* `bbknn` (v1.6.0): Large-scale batch correction.
* `omicverse` (v1.6.10): Pearson residuals scaling and cNMF decomposition.
* `infercnvpy` (v0.5.0): CNV-based malignant cell identification.

### **R Environment (v4.3+)**
* `Seurat` (v4.3.0): Core framework.
* `spacexr` (RCTD v2.2.0): Cell type deconvolution.
* `mistyR` (v1.10.0): Multi-view spatial relationship inference.
* `CellChat` (v1.6.1): Ligand-receptor communication modeling.
* `SPATA2` (v3.1.0): Niche-specific histological gradient screening.

---

## 📝 Author Contributions

**F.W.**, **Z.C.**, and **Z.W.** conceived and designed the study. **F.W.** led the development of methodologies and computational software with support from **R.H.** and **H.L.**. Formal data analysis and visualization were performed by **F.W.**, **R.H.**, **H.L.**, **Y.B.**, **C.Y.**, and **G.Z.**. **F.W.** coordinated the investigation and resource collection with significant contributions from **X.W.**, **W.C.**, **Y.L.**, **Y.N.Z.**, **D.L.**, **Y.Q.**, **J.Z.**, **B.G.**, **C.M.**, **H.D.**, **Y.B.Y.**, **Y.L.Z.**, and **T.S.**. Data curation was managed by **F.W.**, **R.H.**, **H.L.**, and **Y.B.**. **F.W.** wrote the original draft of the manuscript, and **F.W.**, **Z.C.**, and **Z.W.** performed critical revisions and intellectual editing. Overall project supervision and funding acquisition were provided by **Z.C.** and **Z.W.**.

---

## 📝 Citation

If you utilize the GRIT-Atlas code or data, please cite our manuscript:

> **Wang F, et al.** (2025). *A Comprehensive Treatment-Induced Resistance Atlas of Glioblastoma Reveals a Fibrotic Niche Shielding the Tumor from Immunotherapy.bioRxiv 2025.12.25.696471; doi: https://doi.org/10.64898/2025.12.25.696471*

**Contact**:
**Fei Wang, MD**
Department of Neurosurgery & Brain and Nerve Research Laboratory, The First Affiliated Hospital of Soochow University
Email: [wangfeineu@163.com](mailto:wangfeineu@163.com)

---
**License**: 
* **Code**: [MIT License](https://opensource.org/licenses/MIT).
* **Data**: [Creative Commons Zero (CC0)](https://creativecommons.org/publicdomain/zero/1.0/).
