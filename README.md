# GRIT-Atlas: A Comprehensive Treatment-Induced Resistance Atlas of Glioblastoma

[![DOI](https://img.shields.io/badge/Dryad-10.5061/dryad.fbg79cp95-blue.svg)](https://doi.org/10.5061/dryad.fbg79cp95)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational pipeline for the **Glioblastoma Resistance Insights from Treatment Atlas (GRIT-Atlas)**. Our study utilizes an integrated single-cell and spatially resolved multi-omic approach to deconstruct the cellular and spatial logic of therapy-induced evolution in IDH-wildtype glioblastoma.

## Scientific Highlights
* **Scale and Scope**: We constructed a high-resolution integrated atlas encompassing **978,065 cells** from **296 samples**, capturing the full spectrum of disease evolution from primary diagnosis to recurrence post-standard-of-care and combinatorial immunotherapy.
* **The Spatial Resistance Triad**: Using spatial transcriptomics across 48 patient sections, we identified a core functional unit composed of **cNMF7 (MES-like)** malignant cells, differentiation-arrested **E-MDSCs**, and Type VI Collagen-secreting **myCAFs**.
* **Niche-Specific Resistance**: This triad specifically colonizes the hypoxic microvascular proliferation (**MVP**) and pseudopalisading necrosis (**PAN**) niches.
* **Mechanistic Insight**: We demonstrate that myCAFs function as stromal architects, constructing a fibrotic scaffold via the **Collagen/Fibronectin-CD44 signaling axis** to physically exclude cytotoxic T cells and sustain malignant plasticity.

---

## Repository Organization

The analysis is modularized into four key computational stages:

| Script Name | Analysis Module | Primary Methodology |
| :--- | :--- | :--- |
| `1_Single_Cell_Preprocessing.ipynb` | Atlas Construction | Standardized QC, technical covariate regression, and BBKNN-based batch integration. |
| `2_Malignant_cNMF_Analysis.ipynb` | Malignant Phenotyping | Consensus Non-negative Matrix Factorization (cNMF) to delineate 11 functional states. |
| `3_Spatial_Landscape_Analysis.R` | Spatial Modeling | RCTD-based cell type deconvolution and MistyR multi-view spatial relationship inference. |
| `4_Niche_Validation.R` | Signaling & Gradients | CellChat-mediated ligand-receptor network inference and histological gradient screening. |

---

## Data Access & File Structure

All processed data objects required to reproduce the findings are hosted on **Dryad**:
🔗 **[DOI: 10.5061/dryad.fbg79cp95](https://doi.org/10.5061/dryad.fbg79cp95)**

### 1. Global Integrated Data
* **`IntegratedData.zip`**: The master Seurat/AnnData object for the complete GRIT-Atlas (978,065 cells). Includes batch-corrected UMAP coordinates, treatment status (pGBM, rGBM-R+C, rGBM-ICB, rGBM-ICB+Target), and lineage annotations.

### 2. Lineage-Specific Objects (High-Resolution)
* **`Malignant.qs`**: Details the 11 malignant functional states (cNMF), including the resistance-associated **cNMF7 (MES-like)** and NEU states.
* **`Myeloid.zip`**: High-resolution subclustering of 349,583 myeloid cells (23 subpopulations), highlighting **E-MDSC.ADM.HIF1A** and Mono.CD14.FCGR3A.
* **`Lymphocyte.zip`**: Detailed annotations for 36 lymphoid subsets, including cytotoxic T cells, exhausted T cells (**Tex.03.ALDOA.MIF**), and Tregs.
* **`Stroma.qs`**: Characterization of CAFs, pericytes, and SMCs, specifically the **myCAF.02.COL6** subset.
* **`Endothelial.qs`**: High-resolution mapping of the EC compartment, including immune-resistant subsets (**Vein.VCAM1** and **Arterial.01.DKK2**).

### 3. Spatial Transcriptomics (10x Visium)
* **`GRIT-Atlas_RCTDres.zip`**: Comprehensive deconvolution results (RCTD weight matrices) and MistyR spatial importance values for 48 patient sections.
* **`Merge_RCTDInput.qs`**: The optimized scRNA-seq reference object used as the input for RCTD deconvolution.

### 4. Cell-Cell Communication
* **`CellChat.zip`**: Processed CellChat objects documenting the signaling networks within the **"Spatial Resistance Triad"**. Includes three specific sub-objects:
    * `Secreted`: Paracrine/autocrine signaling.
    * `ECM`: Collagen-CD44/Integrin axes (the fibrotic shield).
    * `CCC`: Cell-cell contact-dependent signaling.

---

## Software Requirements

The GRIT-Atlas pipeline is implemented in a hybrid **Python (v3.10.12)** and **R (v4.3.2)** environment.

### **Python Libraries**
* `scanpy` (v1.10.1), `bbknn` (v1.6.0), `omicverse` (v1.6.10), `infercnvpy` (v0.5.0).

### **R Packages**
* `Seurat` (v4.3.0), `spacexr` (v2.2.0), `mistyR` (v1.10.0), `CellChat` (v1.6.1), `SPATA2` (v3.1.0).

---

## Citation

If you utilize the GRIT-Atlas code or data, please cite our manuscript:

> **Wang F, et al.** (2025). *A Comprehensive Treatment-Induced Resistance Atlas of Glioblastoma Reveals a Fibrotic Niche Shielding the Tumor from Immunotherapy.*

**Contact**:
**Fei Wang, MD**
Department of Neurosurgery, the First Affiliated Hospital of Soochow University
Email: [wangfeineu@163.com](mailto:wangfeineu@163.com)

---
**License**: Code (MIT) | Data (CC0).
