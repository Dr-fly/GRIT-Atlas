################################################################################
# Module 4: Intercellular Communication and Spatial Niche Validation
# Focus: The Collagen/Fibronectin-CD44 axis within the Spatial Resistance Triad
# (cNMF7 Malignant Cells + E-MDSCs + myCAF.02.COL6)
################################################################################

library(CellChat)
library(Seurat)
library(tidyverse)
library(qs)
library(patchwork)
library(ComplexHeatmap)

# ==============================================================================
# 1. Data Integration and Preprocessing for CellChat
# ==============================================================================

# Load sub-lineage objects (Malignant, Stroma, and Myeloid)
# Focus specifically on the 'Resistance Triad' subpopulations
scobj_CAF <- qread("Stroma.qs") %>% AdjustGroupLabels()
scobj_CAF <- subset(scobj_CAF, celltype_highres %in% c("myCAF.02.COL6", "myCAF.03.COL1"))

scobj_Mal <- qread("Malignant.qs") %>% AdjustGroupLabels()
scobj_Mal$celltypeForCellchat <- paste0("Malignant_", scobj_Mal$cNMF_cluster)
scobj_Mal <- subset(scobj_Mal, celltypeForCellchat == "Malignant_cNMF_7")

scobj_BMDM <- qread("Myeloid.qs") %>% AdjustGroupLabels()
scobj_BMDM <- subset(scobj_BMDM, celltype %in% c("E_MDSCs.ADM.HIF"))

# Merge triad components into a single object for communication inference
scobj_Triad <- merge(scobj_Mal, y = c(scobj_CAF, scobj_BMDM))
scobj_Triad <- subset(scobj_Triad, group != "Recurrent GBM Immunotherapy") # Focus on SOC and Combo-therapy
scobj_Triad <- NormalizeData(scobj_Triad)

# ==============================================================================
# 2. CellChat Analysis: Inferring the Signaling Infrastructure
# ==============================================================================

# Initialize CellChat object using the curated triad subpopulations
cellchat <- createCellChat(object = scobj_Triad, group.by = "celltypeForCellchat", assay = "RNA")
CellChatDB <- CellChatDB.human

# Analyze specific signaling categories to dissect physical vs. paracrine interactions
# Categories: ECM-Receptor, Secreted Signaling, and Cell-Cell Contact
db_list <- list(
  Total = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
  Secreted = "Secreted Signaling",
  ECM = "ECM-Receptor"
)

run_custom_cellchat <- function(cc_obj, search_tag) {
  cc_obj@DB <- subsetDB(CellChatDB, search = search_tag, key = "annotation")
  cc_obj <- subsetData(cc_obj) %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")
  return(cc_obj)
}

cellchat_ECM <- run_custom_cellchat(cellchat, db_list$ECM)
cellchat_Secreted <- run_custom_cellchat(cellchat, db_list$Secreted)

# ==============================================================================
# 3. Visualization: Signaling Hubs and the Collagen-CD44 Axis
# ==============================================================================

# 

# Plot 1: Fibrotic Scaffold (CAF to Malignant/Myeloid via ECM)
pdf("Results/CellChat_Triad_ECM_Bubble.pdf", height = 8, width = 4)
netVisual_bubble(cellchat_ECM, 
                 sources.use = c("myCAF.02.COL6", "myCAF.03.COL1"), 
                 targets.use = c("Malignant_cNMF_7", "E_MDSCs.ADM.HIF"), 
                 remove.isolate = FALSE,
                 title.name = "myCAF-driven Pro-tumorigenic ECM Signaling")
dev.off()

# Plot 2: Paracrine Feedback (Malignant/Myeloid to CAF via Secreted factors)
pdf("Results/CellChat_Triad_Secreted_Bubble.pdf", height = 6, width = 4)
netVisual_bubble(cellchat_Secreted, 
                 sources.use = c("Malignant_cNMF_7", "E_MDSCs.ADM.HIF"), 
                 targets.use = c("myCAF.02.COL6", "myCAF.03.COL1"), 
                 remove.isolate = FALSE,
                 title.name = "Reciprocal Feedback to Stroma")
dev.off()
