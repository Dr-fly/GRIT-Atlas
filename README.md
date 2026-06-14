# GRIT-Atlas: A Single-Cell and Spatial Atlas of MVP-PAN Evolution Reveals a Desmoplastic Dependency to Overcome Glioblastoma Resistance

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18009414.svg)](https://zenodo.org/records/18009414)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational pipeline and resource links for our study, **"A Single-Cell and Spatial Atlas of MVP-PAN Evolution Reveals a Desmoplastic Dependency to Overcome Glioblastoma Resistance."** By integrating a massive single-cell compendium (GRIT-Atlas) with high-plex spatial molecular imaging (CosMx SMI), we deconstruct the spatiotemporal evolution of therapy-induced resistance in IDH-wildtype glioblastoma.

## 🌟 Scientific Highlights
* **Scale and Scope**: We constructed the comprehensive GRIT-Atlas encompassing **978,065 single cells** from 296 samples, seamlessly integrated with **high-plex CosMx SMI spatial transcriptomics across 406,689 cells** to achieve true single-cell spatial resolution.
* **Deconstructing MVP-PAN Evolution**: We unmasked a deterministic spatiotemporal continuum proving that Pseudopalisading Necrosis (PAN) is the direct pathophysiological consequence of dysfunctional Microvascular Proliferation (MVP).
* **The COL6A1-CD44 Structural Sanctuary**: Stromal cancer-associated fibroblasts (CAFs) act as the primary engine pumping out global extracellular COL6A1, creating a mechanical matrix scar that engages CD44 on MES-like malignant cells to sustain tumor stemness even during vascular collapse.
* **Preclinical Translation via Repurposing**: We identified **Lacidipine**, an FDA-approved blood-brain barrier-penetrant antihypertensive agent, as a potent inhibitor of CAF-mediated COL6A1 secretion, successfully dismantling the desmoplastic niche and overcoming resistance in vivo.

---

## 📂 Repository Organization

The computational workflow is modularized into five core analytical stages:

| Script Name | Analysis Module | Primary Methodology | 
| :--- | :--- | :--- | 
| `1_Single_Cell_Preprocessing_and_Batch_Correction.ipynb` | GRIT-Atlas Construction | QC, Covariate Regression, BBKNN | 
| `2_Malignant_Functional_State_Characterization_via_cNMF.ipynb` | Malignant Phenotyping | Consensus NMF (cNMF) | 
| `3_Spatial_Landscape_Analysis_and_Niche_Characterization.R` | Visium Spatial Modeling | RCTD, MistyR, and SPATA2 |
| `4_Cell_Communication_and_Spatial_Niche_Validation.R` | Cell-Cell Interactions | CellChat Modeling |
| `5_Spatial_Niche_Architecture_Analysis_For_MESlike_Cells.R` | High-Plex CosMx Spatial Niche Mapping | KNN-Neighborhood Extraciton, NMF, Consensus Clustering | 

---

## 📊 Data Access (Zenodo)

All processed data objects required to reproduce the findings are hosted on **Zenodo**:
🔗 **[Direct Link: https://zenodo.org/records/18009414](https://zenodo.org/records/18009414)**

### 1. Global Single-Cell Data (GRIT-Atlas)
* **`IntegratedData.zip`**: The master Seurat/AnnData object (978,065 cells) including batch-corrected UMAP coordinates, treatment metadata, and lineage annotations.
* **`Sub-lineage Objects`**: High-resolution specialized `.qs` and `.zip` objects for the Malignant (including cNMF states), Myeloid, Lymphoid, Stromal, and Endothelial compartments.

### 2. Visium Spatial & Communication Context
* **`GRIT-Atlas_RCTDres.zip`**: RCTD weight matrices and MistyR spatial importance values for 48 glioblastoma patient Visium sections.
* **`CellChat.zip`**: Processed interactome objects validating the **COL6A1-CD44 signaling axis**.

### 3. High-Plex Spatial Molecular Imaging (CosMx SMI)
* **`Annotated_SeuratV5.zip`**: Sample-specific Seurat v5 objects for the high-plex spatial transcriptomics data, fully annotated with global cell-type identities.
* **`Annotated_SeuratV5_MESlike_NicheAnnotation.zip`**: Sample-specific Seurat v5 objects focused exclusively on the MES-like malignant domains, deeply annotated with the identified spatial architectural niches (S1–S3).

### 4. In Vivo Validation scRNA-seq
* **`GL261_scRNAseq.zip`**: The standard expression files for the in vivo orthotopic mouse glioblastoma cohort evaluating the therapeutic efficacy of the Lacidipine triple-combination regimen.

---

## 💻 Software Requirements

### **Python Environment (v3.10+)**
* `scanpy` (v1.10.1): QC and normalization.
* `bbknn` (v1.6.0): Large-scale batch correction.
* `omicverse` (v1.6.10): Pearson residuals scaling and cNMF.
* `infercnvpy` (v0.5.0): CNV-based malignant cell inference.

### **R Environment (v4.3+)**
* `Seurat` (v5.0.0+): Core framework (Essential for CosMx SMI data).
* `spacexr` (v2.2.0): Cell type deconvolution (RCTD).
* `mistyR` (v1.10.0): Multi-view spatial relationship inference.
* `CellChat` (v1.6.1): Ligand-receptor communication modeling.
* `RcppML`(v0.3.7) & `FNN`(v1.1.4.1): Optimized NMF and K-Nearest Neighbors spatial niche extraction.

---

## 👥 Author Contributions (CRediT)

* **Conceptualization:** Fei Wang, Zhouqing Chen, Zhong Wang.
* **Data curation:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang.
* **Formal analysis:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang, Chen Yang, Yuhan Bai.
* **Funding acquisition:** Zhouqing Chen, Zhong Wang, Fei Wang.
* **Investigation:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang, Chen Yang, Yuhan Bai, Wenqian Cao, Yue Lu, Guangling Xu, Haohao Qiu, Hongyi Ling, Dengfeng Lu, Youjia Qiu, Juyi Zhang, Bixi Gao.
* **Methodology:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang, Yanbo Yang, Ting Sun.
* **Project administration:** Zhouqing Chen, Zhong Wang.
* **Resources:** Zhouqing Chen, Zhong Wang, Ting Sun.
* **Software:** Xin Wu, Run Huang, Chen Yang.
* **Supervision:** Zhouqing Chen, Zhong Wang.
* **Validation:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang, Wenqian Cao, Yue Lu.
* **Visualization:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang.
* **Writing – original draft:** Fei Wang, Xin Wu, Guozheng Zhao, Run Huang.
* **Writing – review & editing:** Zhouqing Chen, Zhong Wang, Fei Wang, Yanbo Yang, Ting Sun.

---

## 📝 Citation

If you utilize the code or data from this repository, please cite our manuscript. Note: A preliminary version of this work highlighting the initial single-cell atlas construction has been deposited as a preprint:

> **Wang F, et al.** (2025). *A Comprehensive Treatment-Induced Resistance Atlas of Glioblastoma Reveals a Fibrotic Niche Shielding the Tumor from Immunotherapy.* bioRxiv. doi: [10.64898/2025.12.25.696471](https://doi.org/10.64898/2025.12.25.696471)

**Contact**:  
**Fei Wang, MD** Department of Neurosurgery & Brain and Nerve Research Laboratory, The First Affiliated Hospital of Soochow University  
Email: [wangfeineu@163.com](mailto:wangfeineu@163.com)  

---
**License**:  
* **Code**: [MIT License](https://opensource.org/licenses/MIT).
* **Data**: [Creative Commons Zero (CC0)](https://creativecommons.org/publicdomain/zero/1.0/).
