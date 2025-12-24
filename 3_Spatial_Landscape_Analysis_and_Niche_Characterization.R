################################################################################
# GRIT-Atlas: Complete Spatial Analysis Pipeline
# Methodology: RCTD Deconvolution & MistyR Spatial Modeling
# Focus: Identification of the Spatial Resistance Triad (cNMF7 + E-MDSC + myCAF)
################################################################################

# Load essential libraries
library(spacexr)      
library(mistyR)       
library(future)     
library(Seurat)      
library(tidyverse)    
library(distances)  
library(ComplexHeatmap) 
library(circlize)      
library(pheatmap)    
library(RColorBrewer) 

# Configure global parallel processing strategy
# Allocated 96 workers to handle the intensive computational load of MistyR importance sampling
future::plan(future::multisession(workers = 96))

################################################################################
# PART 1: Dataset and Sample Definitions (48 Patient Samples)
################################################################################

# Define discovery and validation cohorts from 10x Visium platforms
DatasetIDs <- c("CancerCell_Ravi_2023", "NatCancer_Mei_2023", "NatCom_Ren_2023", 
                "Science_Cristina_2025", "Cell_Greenwald_2024")

# Cohort-specific patient identifiers
Ravi_Samples <- c("#UKF242_T_ST", "#UKF243_T_ST", "#UKF248_T_ST", "#UKF251_T_ST", 
                  "#UKF255_T_ST", "#UKF259_T_ST", "#UKF260_T_ST", "#UKF262_T_ST", 
                  "#UKF265_T_ST", "#UKF266_T_ST", "#UKF269_T_ST", "#UKF275_T_ST", 
                  "#UKF296_T_ST", "#UKF304_T_ST", "#UKF313_T_ST", "#UKF334_T_ST")

Mei_Samples <- c("Pt1", "Pt2", "Pt3", "Pt4", "Pt5", "Pt6", "Pt9", "Pt10", 
                 "Pt11", "Pt12", "Pt13", "Pt14", "Pt15", "Pt16", "Pt17", "Pt18")

Ren_Samples <- c("GBM2", "GBM3", "GBM5_1")

Cristina_Samples <- c("GBM03422", "GBM25526")

Greenwald_Samples <- c("mgh258", "zh1007inf", "zh1007nec", "zh1019inf", "zh1019t1", 
                       "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1", 
                       "zh916bulk", "zh916inf", "zh916t1")

################################################################################
# PART 2: Spatial Deconvolution using RCTD
################################################################################

# Function to estimate cell type proportions using the GRIT-Atlas reference
RunRCTD <- function(DatasetID, PatientID) {
  # Load the GRIT-Atlas single-cell reference
  reference <- qs::qread("GRIT_Atlas/Spatial/scRCTD/Merge_RCTDInput.qs")
  
  # Load 10x Visium spatial query data
  load(paste("GRIT_Atlas/Spatial/", DatasetID, "/RCTDInput/", PatientID, "_RCTDInput.Rdata", sep=""))
  
  # Execute RCTD in 'doublet' mode to resolve cell type mixtures within spatial spots
  RCTD <- create.RCTD(query, reference, max_cores = 96)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Synchronize spatial spots and integrate weights into Seurat metadata
  seu_vis <- AddMetaData(seu_vis, metadata = RCTD@results$results_df)
  weights.df <- as.data.frame(as.matrix(RCTD@results$weights))
  commonspots <- intersect(colnames(seu_vis), rownames(weights.df))
  seu_vis <- subset(seu_vis, cells = commonspots)
  seu_vis <- AddMetaData(seu_vis, metadata = weights.df)
  
  # Save processed objects
  save(seu_vis, RCTD, 
       file = paste("GRIT_Atlas/Spatial/", DatasetID, "/", PatientID, "_RCTDres.Rdata", sep=""))
}

# Batch processing for deconvolution
for (ID in Ravi_Samples)      RunRCTD("CancerCell_Ravi_2023", ID)
for (ID in Mei_Samples)       RunRCTD("NatCancer_Mei_2023", ID)
for (ID in Ren_Samples)       RunRCTD("NatCom_Ren_2023", ID)
for (ID in Cristina_Samples)  RunRCTD("Science_Cristina_2025", ID)

################################################################################
# PART 3: Inter-cellular Relationship Inference (MistyR)
################################################################################

# Function to quantify importance of cNMF7 as a predictor of TME composition
RunMistyR <- function(DatasetID, PatientID) {
  # Load deconvolution results
  load(paste("GRIT_Atlas/Spatial/", DatasetID, "/", PatientID, "_RCTDres.Rdata", sep=""))
  
  # Extract composition matrix (Target variables)
  metadata <- seu_vis@meta.data
  composition <- as_tibble(metadata[, 13:ncol(metadata)])
  rownames(composition) <- colnames(seu_vis)
  
  # Calculate characteristic neighborhood radius (Paraview radius)
  geometry <- GetTissueCoordinates(seu_vis, cols = c("imagerow", "imagecol"), scale = NULL)
  geom_dist <- as.matrix(distances(geometry))  
  dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
  paraview_radius <- ceiling(mean(dist_nn + sd(dist_nn)))
  
  # Construct Intraview (local) and Paraview (neighborhood) views
  misty_views <- create_initial_view(composition) %>%
    add_paraview(geometry, l = paraview_radius, family = "gaussian")
  
  # Run MISTy and aggregate results
  output_folder <- paste0("result/", PatientID, "_misty_pipeline")
  run_misty(misty_views, output_folder)
  misty_results <- collect_results(output_folder)
  
  # Save spatial importance values
  save(misty_results, file = paste("GRIT_Atlas/Spatial/RCTDres/", 
                                   DatasetID, "/", PatientID, "_MistyRres.Rdata", sep=""))
  unlink(output_folder, recursive = TRUE)
}

# Batch processing for spatial modeling
for (ID in Ravi_Samples)      RunMistyR("CancerCell_Ravi_2023", ID)
for (ID in Mei_Samples)       RunMistyR("NatCancer_Mei_2023", ID)
for (ID in Ren_Samples)       RunMistyR("NatCom_Ren_2023", ID)
for (ID in Cristina_Samples)  RunMistyR("Science_Cristina_2025", ID)
for (ID in Greenwald_Samples) RunMistyR("Cell_Greenwald_2024", ID)

################################################################################
# PART 4: Quantitative Meta-Analysis & Visualization
################################################################################

# Aggregate results focusing on the cNMF7 malignant subpopulation
locations <- c(paste("CancerCell_Ravi_2023", Ravi_Samples, sep="/"),
               paste("NatCancer_Mei_2023", Mei_Samples, sep="/"),
               paste("NatCom_Ren_2023", Ren_Samples, sep="/"),
               paste("Science_Cristina_2025", Cristina_Samples, sep="/"),
               paste("Cell_Greenwald_2024", Greenwald_Samples, sep="/"))

temp0 <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(temp0) <- c("view", "Predictor", "Target", "Importance", "nsamples", "DatasetID", "PatientID")

for (location in locations) {
  load(paste("GRIT_Atlas/Spatial/RCTDres/", location, "_MistyRres.Rdata", sep=""))
  df <- misty_results[["contributions.stats"]]
  temp_stats <- df[order(df$fraction, decreasing = T),] %>% head(93)
  
  temp_imp <- misty_results$importances.aggregated %>% filter(Predictor == "Malignant_cNMF_7")
  temp_imp$tool <- paste(temp_imp$Target, temp_imp$view, sep=" ")
  
  temp_final <- temp_imp[temp_imp$tool %in% temp_stats$tool,] %>% arrange(-Importance)
  temp_final$ID <- paste(strsplit(location, split = "/")[[1]][1], strsplit(location, split = "/")[[1]][2], sep="-")
  temp0 <- rbind(temp0, temp_final %>% filter(Target != "Malignant_cNMF_7"))
}

# Categorize and pivot for Niche Characterization
temp <- as.data.frame(temp0)
temp$group <- ifelse(temp$PatientID %in% c("Pt2","Pt11","Pt15","Pt16","Pt17","Pt18"), "rGBM(ICB+T)",
                     ifelse(temp$PatientID %in% c("Pt3","Pt6","Pt12"), "rGBM(C+R)", "pGBM"))
temp$response <- ifelse(temp$PatientID %in% c("Pt2","Pt11","Pt16","Pt17"), "NR",
                        ifelse(temp$PatientID %in% c("Pt15","Pt18"), "R", "Unknown"))

wide_data <- temp %>%
  dplyr::select(Target, ID, Importance) %>%
  tidyr::pivot_wider(names_from = ID, values_from = Importance, values_fill = NA) %>%
  tibble::column_to_rownames(var = "Target")

# Define target features
features_to_plot <- c(
  "Malignant_cNMF_1", "Stroma_apCAF_CD74", "E_MDSCs_ADM_HIF", "Malignant_cNMF_5",
  "Stroma_PC_01_RGS5", "EC_Cap_02_CA4", "M_MDSCs_CXCL1_CXCL2_CXCL3_FCN1_IL10", "TAM_03_RPL",
  "Mono_CD14_FCGR3A", "B_Bm_02_HSPA1A", "Stroma_prolCAF_TOP2A", "EC_Arterial_03_APLN",
  "Malignant_cNMF_4", "TAM_02_APOE", "Malignant_cNMF_6", "Stroma_myCAF_03_COL1",
  "cDC1_02_IRF8_CLEC9A", "B_Bm_01_CD79A", "EC_Cap_04_MKI67", "Malignant_cNMF_8",
  "MG_04_TXNIP_S100A", "Stroma_myCAF_02_COL6", "Stroma_SMC_ACTA2_MYH11", "Malignant_cNMF_3",
  "Mast", "EC_Arterial_02_CXCL12"
)
wide_data <- wide_data[features_to_plot,]

# Binarization and frequency analysis (10% to 50% thresholds)
GetSpecificRatio <- function(data, percentile) {
  convert_to_binary <- function(input_data, p) {
    binary_data <- input_data
    for (i in 1:ncol(binary_data)) {
      col_values <- as.numeric(binary_data[, i])
      n_top <- round(length(col_values) * p)
      binary_data[, i] <- 0
      binary_data[order(col_values, decreasing = TRUE)[1:n_top], i] <- 1
    }
    return(binary_data)
  }
  binary_mat <- convert_to_binary(data, percentile)
  count_ones <- rowSums(binary_mat == 1, na.rm = TRUE)
  data.frame(CellType = rownames(binary_mat), Count_of_1 = count_ones,
             Percentage = round(count_ones / ncol(binary_mat) * 100, 1)) %>% arrange(-Count_of_1)
}

# Final Frequency Heatmap Generation
freq_list <- lapply(seq(0.1, 0.5, 0.1), function(p) GetSpecificRatio(wide_data, p))
percent_matrix <- do.call(cbind, lapply(freq_list, function(x) x$Percentage))
counts_matrix <- do.call(cbind, lapply(freq_list, function(x) x$Count_of_1))
colnames(percent_matrix) <- colnames(counts_matrix) <- seq(0.1, 0.5, 0.1)
rownames(percent_matrix) <- rownames(counts_matrix) <- freq_list[[1]]$CellType

pdf("GRIT_Atlas/Results/Figure/Figure6/FrequencyHM.pdf", width = 5, height = 15)
pheatmap(percent_matrix, display_numbers = counts_matrix, cluster_rows = T, cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Recurrent cNMF7 Spatial Associations Across Cohorts")
dev.off()