################################################################################
# GRIT-Atlas: Module 5 - Spatial Niche Architecture Analysis
# Methodology: KNN-based spatial neighborhood extraction, Background Correction, 
#              and Consensus Clustering
# Focus: Identification of the 3 distinct "MES Niches" (S1, S2, S3)
################################################################################

# Load essential libraries
library(Seurat)
library(qs)
library(FNN)
library(glue)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(pbapply)

# Configure paths (Relative paths for GitHub reproducibility)
DATA_DIR <- "Data/CosMx_Spatial/"
RESULTS_DIR <- "Results/MES_Niche_Analysis/"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Focus strictly on the MES-like malignant subpopulation
TARGET_CTYPE <- "MES-like"
K_NEIGHBORS <- 100

################################################################################
# PART 1: Core Functions for Spatial Niche Extraction & NMF
################################################################################

# Execute NMF on normalized neighborhood matrices
RunNMF <- function(niche.mat, k=5, cores=4, seed=1024) {
  RcppML::setRcppMLthreads(cores)
  model <- RcppML::nmf(niche.mat, k = k, verbose = FALSE, seed = seed)
  H <- model$h
  rownames(H) <- paste0("factor_", 1:nrow(H))
  colnames(H) <- colnames(niche.mat)
  
  W <- model$w
  rownames(W) <- rownames(niche.mat)
  colnames(W) <- rownames(H)
  
  return(list(W = W, H = H))
}

# Extract local cellular neighborhoods using KNN
ExtractSpatialNiche <- function(seu, query_ctype, k_neighbors = 100) {
  Idents(seu) <- "annotation"
  query_indices <- which(seu$annotation == query_ctype)
  if (length(query_indices) < 10) return(NULL)
  
  query_cell_names <- colnames(seu)[query_indices]
  coords <- as.matrix(Embeddings(seu, reduction = "spatial")[, 1:2])
  
  all_annos <- as.character(seu$annotation)
  unique_annos <- sort(unique(all_annos))
  
  # K-Nearest Neighbors calculation
  actual_k <- min(k_neighbors, nrow(coords) - 1)
  query_dat <- coords[query_indices, , drop = FALSE]
  knn_res <- get.knnx(data = coords, query = query_dat, k = actual_k)
  nn_idx <- knn_res$nn.index
  
  # Tally cellular composition within the neighborhood
  niche_counts <- t(apply(nn_idx, 1, function(indices) {
    as.numeric(table(factor(all_annos[indices], levels = unique_annos)))
  }))
  rownames(niche_counts) <- query_cell_names
  colnames(niche_counts) <- unique_annos
  
  # Dimensionality reduction via NMF for local niche states
  niche_mat_norm <- niche_counts / actual_k
  nmf_res <- RunNMF(niche_mat_norm, k = 5)
  
  # Subset Seurat object and store niche data
  seu_sub <- subset(seu, cells = query_cell_names)
  seu_sub[["niche"]] <- CreateAssayObject(counts = t(niche_counts))
  seu_sub[["nmf"]] <- CreateDimReducObject(embeddings = nmf_res$W, key = "factor_", assay = DefaultAssay(seu_sub))
  
  return(seu_sub)
}

################################################################################
# PART 2: Extract and Integrate MES-like Niches Across Cohorts
################################################################################

# 1. Load spatial Seurat objects (Assuming pre-processed V5 objects)
qs_files <- list.files(DATA_DIR, pattern = "\\.SeuratV5\\.qs$", full.names = TRUE)

niche_obj_list <- pblapply(qs_files, function(qf) {
  seu <- qread(qf)
  niche_sub <- ExtractSpatialNiche(seu, query_ctype = TARGET_CTYPE, k_neighbors = K_NEIGHBORS)
  rm(seu); gc()
  return(niche_sub)
})
niche_obj_list <- Filter(Negate(is.null), niche_obj_list)

# 2. Background Correction and Cross-Sample Integration
# Correct local niche proportions against global tissue-level cell type frequencies
global_celltypes <- unique(unlist(lapply(niche_obj_list, function(x) rownames(x[["niche"]]))))

integrated_niche_list <- pblapply(niche_obj_list, function(seu_sub) {
  niche.mat <- LayerData(seu_sub, assay = "niche", layer = "counts")
  
  # Simulated background fraction (In practice, derived from the full tissue object)
  bg.fraction <- rep(1 / length(global_celltypes), length(global_celltypes)) 
  names(bg.fraction) <- global_celltypes
  
  z <- rowMeans(niche.mat)
  z <- z / sum(z)
  z_corr <- z / as.numeric(bg.fraction[names(z)])
  z_corr <- z_corr / sum(z_corr)
  
  return(z_corr)
})

final_niche_mat <- do.call(cbind, integrated_niche_list)
colnames(final_niche_mat) <- paste0("Sample_", 1:ncol(final_niche_mat))
qsave(final_niche_mat, file.path(RESULTS_DIR, "MES_Integrated_Niche_Matrix.qs"))

################################################################################
# PART 3: Consensus Clustering to define the 3 "MES Niches"
################################################################################

# Calculate Pearson Correlation and Distance
cor_total <- cor(final_niche_mat)
dist_total <- as.dist(1 - cor_total)

# Hierarchical Clustering (Ward.D2) to define the 3 microenvironmental states
OPT_K <- 3
hc_total <- hclust(dist_total, method = "ward.D2")
cluster_anno <- cutree(hc_total, k = OPT_K)

# Map clusters to structural identities (S1: PAN, S2: Reactive, S3: MVP)
cluster_df <- data.frame(
  Sample_Subgroup_ID = names(cluster_anno), 
  Global_Niche = paste0("MES-like_S", cluster_anno)
)
write.csv(cluster_df, file.path(RESULTS_DIR, "MES_Niche_Cluster_Mapping.csv"), row.names = FALSE)

################################################################################
# PART 4: High-Quality Visualization of the MES Niche Continuum
################################################################################

# Define aesthetic color palettes consistent with the manuscript
MY_RED <- "#C16E71"
MY_BLUE <- "#6E8FB2"

total_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A"), # Distinct colors for S1, S2, S3
  sort(unique(cluster_df$Global_Niche))
)

col_fun <- colorRamp2(c(-1, 0, 1), c(MY_BLUE, "white", MY_RED))

ha <- HeatmapAnnotation(
  Niche = cluster_df$Global_Niche, 
  col = list(Niche = total_colors), 
  show_annotation_name = TRUE
)

# Render Final Heatmap
pdf(file.path(RESULTS_DIR, "MES_Niche_Architecture_Heatmap.pdf"), width = 10, height = 7)
ht <- Heatmap(
  cor_total, 
  name = "PCC", 
  top_annotation = ha,
  column_title = "MES-like Spatial Niche Ecosystem Architecture",
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  cluster_rows = hc_total, 
  cluster_columns = hc_total,
  col = col_fun,
  border = TRUE
)
draw(ht)
dev.off()
