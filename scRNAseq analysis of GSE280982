
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("patchwork")

# Install Bioconductor packages for advanced integration
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "batchelor",  # FastMNN integration (used by Seurat's IntegrateLayers)
  "scran"       # Additional normalization methods
))

# Install Harmony for integration
install.packages("harmony")
install.packages("rlang")
# Install additional visualization and utility packages
install.packages(c(
  "remotes",         # Install R packages stored in GitHub
  "ggrepel",         # Non-overlapping text labels
  "RColorBrewer",    # Color palettes
  "viridis"          # Perceptually uniform color scales
))

# Wrapper functions for integration methods
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library("remotes")
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
remotes::install_github('satijalab/seurat-wrappers')
# Install future for parallel processing (used by RPCA integration)
install.packages("future")

# Install FNN for k-nearest neighbor calculations (used in quality metrics)
install.packages("FNN")

# Install cluster for clustering quality metrics (silhouette scores)
install.packages("cluster")

# Install reshape2 for data reshaping (used in cluster quality assessment)
install.packages("reshape2")
install.packages("Seurat")

###///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
##Core single-cell analysis
library(SeuratObject)

# Integration methods
library(harmony)
library(batchelor)
library(SeuratWrappers)

# Visualization and data manipulation
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Quality metrics and utilities
library(FNN)              # K-nearest neighbor calculations for mixing metrics
library(cluster)          # Silhouette scores for clustering quality
library(reshape2)         # Data reshaping for visualization

# Parallel processing
library(future)
options(future.globals.maxSize = 20 * 1024^3)  # Increase to 20GB for large datasets

library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(glmGamPoi)
library(ggpubr)
library(scales)
##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ==============================================================================
# STEP 1: Define the "Bulletproof" Loading Function
# (This wraps the manual code that worked for you)
# ==============================================================================
load_geo_sample <- function(folder_path, prefix, project_name) {
  
  message(paste0("--- Processing: ", project_name, " ---"))
  
  # Construct filenames automatically based on the prefix
  # Note: Adjust pattern if your files are named differently (e.g., "genes.tsv" instead of "features.tsv")
  mtx_path  <- file.path(folder_path, paste0(prefix, "matrix.mtx.gz"))
  feat_path <- file.path(folder_path, paste0(prefix, "features.tsv.gz"))
  bc_path   <- file.path(folder_path, paste0(prefix, "barcodes.tsv.gz"))
  
  # Check if files exist
  if (!file.exists(mtx_path) || !file.exists(feat_path) || !file.exists(bc_path)) {
    stop(paste("Missing files for prefix:", prefix, "\nCheck your spelling!"))
  }
  
  # Manual Load (The method that worked for you)
  mat <- readMM(gzfile(mtx_path))
  features <- read.delim(gzfile(feat_path), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(bc_path), header = FALSE, stringsAsFactors = FALSE)
  
  # Assign Names
  # Handle duplicate gene names
  rownames(mat) <- make.unique(features$V2) 
  colnames(mat) <- barcodes$V1
  
  # Create Seurat Object
  obj <- CreateSeuratObject(counts = mat, project = project_name)
  return(obj)
}

# ==============================================================================
# STEP 2: Define Your Samples
# ==============================================================================

# Set your folder path
data_dir <- "/GSE280982_RAW"

# UPDATE THIS LIST with the exact filenames from your folder!
# format: "Filename_Prefix_" = "Short_Name_You_Want"
samples_info <- list(
  # Example: "GSM8608318_HyPR-HN_01_1_" is the file prefix, "Control_1" is the label
  "GSM8608318_HyPR-HN_01_1_" = "prerad_1",
  "GSM8608319_HyPR-HN_01_2_" = "rad_1",  # Check your actual filename for this!
  "GSM8608320_HyPR-HN_01_3_" = "rad6w_1",
  "GSM8608323_HyPR-HN_02_1_" = "prerad_2",      # Check your actual filename!
  "GSM8608328_HyPR-HN_02_2_" = "rad_2",
  "GSM8608330_HyPR-HN_02_3_" = "rad6w_2",
  "GSM8608338_HyPR-HN_03_1_" = "prerad_3",
  "GSM8608340_HyPR-HN_03_2_" = "rad3",
  "GSM8608344_HyPR-HN_04_1_" = "prerad4",
  "GSM8608346_HyPR-HN_04_2_" = "rad4",
  "GSM8608348_HyPR-HN_04_3_" = "rad6w_4"
  
  
)

# ==============================================================================
# STEP 3: Load All Samples
# ==============================================================================
seurat_list <- list()

for (prefix in names(samples_info)) {
  sample_name <- samples_info[[prefix]]
  
  # Run the function
  seurat_list[[sample_name]] <- load_geo_sample(data_dir, prefix, sample_name)
}

# ==============================================================================
# STEP 4: Merge into One Object
# ==============================================================================
message("Merging all samples...")

# We take the first object and merge the rest into it
combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list), # This adds "Control_1_" to barcode names so they don't clash
  project = "HNSCC_Combined"
)

# ==============================================================================
# STEP 5: Verify
# ==============================================================================
print(combined_seurat)
table(combined_seurat$orig.ident) # Shows how many cells are in each sample

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////QC
# 1. Calculate Mitochondrial Percentage for all cells
# (This looks for genes starting with "MT-")
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")

# 2. Visualize QC metrics
# We group by "orig.ident" so you can compare the 4 samples side-by-side
VlnPlot(combined_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "orig.ident", 
        ncol = 3,
        pt.size = 0) # pt.size=0 hides the dots so the plot is cleaner


# 1. Extract the metadata
metadata <- combined_seurat@meta.data

# 2. Calculate averages and medians for each sample
qc_stats <- metadata %>%
  group_by(orig.ident) %>%
  summarise(
    # Count how many cells are in each group
    n_Cells = n(),
    
    # Average and Median Genes detected (nFeature)
    Mean_Genes = mean(nFeature_RNA),
    Median_Genes = median(nFeature_RNA),
    
    # Average and Median Total Molecules (nCount)
    Mean_Counts = mean(nCount_RNA),
    Median_Counts = median(nCount_RNA),
    
    # Average and Median Mitochondrial %
    Mean_MT = mean(percent.mt),
    Median_MT = median(percent.mt)
  )

# 3. Print the table nicely
print(qc_stats)

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#######Filtering cells
###########
##########
# 1. Define your cutoffs (Variables make it easier to change later)
min_feature <- 200
max_feature <- 5000
min_count   <- 1000
max_count   <- 20000
max_mt      <- 8

# 2. Apply the subset command
# The subset() function uses logical AND (&) by default
seurat_filtered <- subset(combined_seurat, 
                          subset = nFeature_RNA > min_feature & 
                            nFeature_RNA < max_feature & 
                            nCount_RNA   > min_count & 
                            nCount_RNA   < max_count & 
                            percent.mt   < max_mt)

# 3. Compare Before vs. After
message("--- Cell Counts ---")
message("Original Cells: ", ncol(combined_seurat))
message("Filtered Cells: ", ncol(seurat_filtered))
message("Cells Removed:  ", ncol(combined_seurat) - ncol(seurat_filtered))

# 4. Check the new breakdown by sample
table(seurat_filtered$orig.ident)

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#####filter gene
counts_matrix <- LayerData(seurat_filtered, layer = "counts")
gene_detection <- rowSums(counts_matrix > 0)

gene_qc <- data.frame(
  gene = rownames(seurat_filtered),
  n_cells_detected = gene_detection,
  pct_cells_detected = (gene_detection / ncol(seurat_filtered)) * 100,
  is_mt = grepl("^MT-", rownames(seurat_filtered)),
  is_ribo = grepl("^RP[SL]", rownames(seurat_filtered)),
  is_hb = grepl("^HB[AB]", rownames(seurat_filtered))
)

cat("Total genes:", nrow(gene_qc), "\n")
cat("MT genes:", sum(gene_qc$is_mt), 
    "| Ribo genes:", sum(gene_qc$is_ribo), 
    "| Hb genes:", sum(gene_qc$is_hb), "\n\n")

# Visualize gene detection
library(ggplot2)
p9 <- ggplot(gene_qc, aes(x = pct_cells_detected)) +
  geom_histogram(bins = 50, fill = "#118AB2", alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(
    title = "Gene Detection Across Cells",
    subtitle = "Red line: 0.1% of cells threshold",
    x = "% of Cells Expressing Gene (log scale)",
    y = "Number of Genes"
  ) +
  theme_classic()

cat("→ EXAMINE plots/06_gene_detection.png\n")
cat("   Choose threshold based on where detection drops off\n\n")

# Set filtering threshold - ADJUST BASED ON YOUR DATA
# Common approaches:
#   - 0.1% of cells (very lenient, keeps most genes)
#   - 1% of cells (moderate, standard approach)
#   - 3 cells minimum (absolute count, conservative)

min_pct_cells <- 0.1  # Gene must be in ≥0.1% of cells
# Alternatively: min_cells <- 3  # Gene must be in ≥3 cells

cat("Setting gene filter threshold:\n")
cat("  Minimum detection: ≥", min_pct_cells, "% of cells\n")
cat("  (Equivalent to ≥", ceiling(ncol(seurat_filtered) * min_pct_cells / 100), 
    "cells with current dataset)\n\n")

cat("Threshold options:\n")
cat("  • 0.1% of cells: Lenient (keeps rare cell type markers)\n")
cat("  • 1% of cells: Standard (balances detection vs noise)\n")
cat("  • 3-10 cells minimum: Conservative (removes very rare genes)\n\n")

# Filter genes
genes_to_keep <- (gene_qc$pct_cells_detected >= min_pct_cells) & !gene_qc$is_hb

cat("Genes passing filter:", sum(genes_to_keep), "/", nrow(gene_qc), "\n")
cat("Genes removed:\n")
cat("  Low detection:", sum(gene_qc$pct_cells_detected < min_pct_cells), "\n")
cat("  Hemoglobin:", sum(gene_qc$is_hb), "\n\n")

seurat_filtered <- seurat_filtered[gene_qc$gene[genes_to_keep], ]

cat("After gene filtering:", nrow(seurat_filtered), "genes remaining\n")

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

####normalize, scale and pca
#--- PART 1: Standard Preprocessing ---
  # We must Normalize, Scale, and Run PCA before we can integrate.

# 1. Ensure you have the percentage calculated
seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_filtered, pattern = "^MT-")

# Run SCTransform
# vars.to.regress = "percent.mt" tells Seurat to ignore mitochondrial noise
seurat_sct <- SCTransform(seurat_filtered, 
                          vars.to.regress = "percent.mt", 
                          verbose = FALSE)
message("Step 2: Finding Variable Features...")
# Identify the top 2000 genes that define biological differences
seurat_sct <- FindVariableFeatures(seurat_sct, selection.method = "vst", nfeatures = 2000)

message("Step 3: Scaling Data...")
# Centers expression to 0 and variance to 1 (crucial for PCA)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

message("Step 4: Running PCA...")
# Reduces dimensions to prepare for UMAP/Harmony
seurat_sct <- RunPCA(seurat_sct, features = VariableFeatures(object = seurat_sct))

###visualize the scaling
top10 <- head(VariableFeatures(seurat_sct), 10)

# Plot variable features with and without labels
plot3 <- ElbowPlot(seurat_sct) ###elbow plot

plot4 <- DimHeatmap(seurat_sct, dims = 1, cells = 500, balanced = TRUE)###PC1 heatmap

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
###intergration

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/integration_comparison", showWarnings = FALSE)
dir.create("plots/clustering", showWarnings = FALSE)
dir.create("integrated_data", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)
# Configure plotting defaults
theme_set(theme_classic(base_size = 12))

current_layers <- Layers(seurat_sct[["RNA"]])
if (length(current_layers) > 1) {
  cat("Layers already split by merge operation (", length(current_layers), " layers)\n", sep = "")
} else {
  # Split layers by sample if not already split
  cat("Splitting layers by sample\n")
  seurat_sct[["RNA"]] <- split(seurat_sct[["RNA"]], f = seurat_sct$sample_id)
}


######integration by harmony
message("Running Harmony...")
seurat_filtered_H <- RunHarmony(seurat_sct, 
                              group.by.vars = "orig.ident", 
                              plot_convergence = FALSE)

#######We run a NEW UMAP based on the corrected 'harmony' reduction
message("Running Harmony UMAP...")
seurat_filtered_Humap <- RunUMAP(seurat_filtered_H, 
                           dims = 1:15, 
                           reduction = "harmony", 
                           reduction.name = "umap.harmony",
                           seed.use = 42)
# We group "prerad_1" & "prerad_2" -> "Prerad"
# We group "Rad_1" & "Rad_2"       -> "Rad"
seurat_filtered_Humap$condition <- case_when(
  grepl("prerad", seurat_filtered_Humap$orig.ident, ignore.case = TRUE) ~ "Prerad",
  grepl("rad6w",  seurat_filtered_Humap$orig.ident, ignore.case = TRUE) ~ "Rad6w",
  grepl("rad",    seurat_filtered_Humap$orig.ident, ignore.case = TRUE) ~ "Rad"
)

# This forces the plots to appear in your preferred chronological order
seurat_filtered_Humap$condition <- factor(seurat_filtered_Humap$condition, 
                                          levels = c("Prerad", "Rad", "Rad6w"))


# Plot 2: Harmony (Watch for mixing)
p2 <- DimPlot(seurat_filtered_Humap, reduction = "umap.harmony", group.by = "orig.ident") + 
  ggtitle("Harmony Integrated") + 
  theme(plot.title = element_text(hjust = 0.5))

###plot side by side naive merge and harmony
p2

###sample split plot of naive merge umap
DimPlot(seurat_filtered_NM, 
        reduction = "umap.naive", 
        split.by = "orig.ident", 
        ncol = 2) +  # Arranges plots in 2 columns
  ggtitle("Naive Integration: Split by Sample")

###sample split plot of harmony umap
DimPlot(seurat_filtered_Humap, 
        reduction = "umap.harmony", 
        split.by = "orig.ident", 
        ncol = 2) + 
  ggtitle("Harmony Integrated: Split by Sample")

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
###test biological preservence
# Function to calculate mixing metric (local inverse Simpson's Index)
# This measures sample diversity in each cell's k-nearest neighborhood

# Verify the grouping worked (You should see counts for just Prerad and Rad)
print("Cell counts per condition:")
print(table(seurat_filtered_Humap$condition))

# --- STEP 2: Define the Distance Calculator ---
calculate_centroid_dist <- function(object, reduction_name) {
  
  # 1. Extract UMAP coordinates
  # We use the reduction you specify (naive or harmony)
  umap_coords <- Embeddings(object, reduction = reduction_name) %>% as.data.frame()
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  # 2. Add the condition info to the coordinates
  umap_coords$condition <- object$condition
  
  # 3. Calculate the "Center of Gravity" (Centroid) for each group
  centers <- umap_coords %>%
    group_by(condition) %>%
    summarise(
      mean_x = mean(UMAP_1),
      mean_y = mean(UMAP_2)
    )
  
  # 4. Calculate Euclidean Distance between Prerad and Rad centers
  # Formula: sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
  x1 <- centers$mean_x[1]; y1 <- centers$mean_y[1] # Prerad Center
  x2 <- centers$mean_x[2]; y2 <- centers$mean_y[2] # Rad Center
  
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}
  

# --- STEP 3: Run and Compare ---

# Calculate for Naive (No correction)
dist_naive <- calculate_centroid_dist(seurat_filtered_NM, "umap.naive")

# Calculate for Harmony (Corrected)
dist_harmony <- calculate_centroid_dist(seurat_filtered_Humap, "umap.harmony")
  
# --- STEP 4: Print Results ---
message("------------------------------------------------")
message("DISTANCE BETWEEN PRERAD AND RAD CENTERS:")
message("------------------------------------------------")
message(paste("Naive Integration (Batch Effect): ", round(dist_naive, 3)))
message(paste("Harmony Integration (Corrected):  ", round(dist_harmony, 3)))

# Calculate improvement
improvement <- round((dist_naive - dist_harmony), 3)
message(paste("Distance Reduced by:              ", improvement))

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////clustering
##findneighbors
seurat_filtered_Humap <- FindNeighbors(seurat_filtered_Humap, 
                                       reduction = "harmony", 
                                       dims = 1:15)

###clustering on harmony intergrated object
# Test multiple resolutions
resolutions <- c(0.1, 0.3, 0.5)

for (res in resolutions) {
  seurat_filtered_Humap <- FindClusters(
    seurat_filtered_Humap,
    resolution = res,
    verbose = FALSE
  )
  
  n_clusters_SCT <- length(unique(seurat_filtered_Humap@meta.data[[paste0("SCT_snn_res.", res)]]))
  cat(sprintf("Resolution %.1f: %d clusters\n", res, n_clusters_SCT))

}
# Rename cluster columns for clarity
colnames(seurat_filtered_Humap@meta.data) <- gsub("RNA_snn_res\\.", "SCT_snn_res.", 
                                             colnames(seurat_filtered_Humap@meta.data))

# Create UMAP plots for each resolution
plot_list <- lapply(resolutions, function(res) {
  cluster_col <- paste0("SCT_snn_res.", res)
  n_clusters <- length(unique(seurat_filtered_Humap@meta.data[[cluster_col]]))
  
  DimPlot(seurat_filtered_Humap, 
          reduction = "umap.harmony",
          group.by = cluster_col,
          label = TRUE,
          label.size = 4,
          pt.size = 0.3) +
    ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
    NoLegend()
})

# Combine plots
combined_resolutions <- wrap_plots(plot_list, ncol = 2)
ggsave("plots/clustering/06_multi_resolution_clustering.png", 
       combined_resolutions, width = 14, height = 14, dpi = 300)
###find clusterts set resolution at 0.1
seurat_filtered_Humap <- FindClusters(seurat_filtered_Humap, resolution = 0.1)

##check the cells number from each sample
table(seurat_filtered_Humap$orig.ident)

###cells number for each group
table(seurat_filtered_Humap$condition)
##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Cluster QC
# For large datasets (>10,000 cells), subsample for silhouette calculation (use 5000 cells in our case)
# Computing distance matrices for all cells is computationally prohibitive
n_cells <- ncol(seurat_filtered_Humap)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  cat("Dataset has", n_cells, "cells - subsampling", max_cells_for_silhouette, "cells for silhouette calculation\n")
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else {
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)
  
  # Get clusters for subsampled cells
  clusters <- as.numeric(seurat_filtered_Humap@meta.data[[cluster_col]][subsample_idx])
  
  # Get integrated coordinates for subsampled cells (use first 30 dimensions)
  coords <- Embeddings(seurat_filtered_Humap, reduction = "umap.harmony")[subsample_idx, ]
  coords <- coords[, 1:min(30, ncol(coords))]
  
  # Calculate silhouette (only if we have at least 2 clusters)
  if (length(unique(clusters)) > 1) {
    dist_matrix <- dist(coords)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  } else {
    NA  # Return NA if only 1 cluster
  }
})

# Create resolution comparison table
resolution_comparison <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(res) {
    cluster_col <- paste0("clusters_res_", res)
    length(unique(seurat_filtered_Humap@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

print(resolution_comparison)  ###looking for silhouette_score at least 0.3

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////visualized integrated cluster treatment group side by side
# 1. Ensure the identities are set to your 0.1 resolution clusters
# (Seurat automatically sets the active identity to the last FindClusters run)
Idents(seurat_filtered_Humap) <- "seurat_clusters"

# 2. Generate the Split Plot
# Colors = Clusters (0, 1, 2...)
# Panels = Conditions (Prerad vs Rad)
DimPlot(seurat_filtered_Humap, 
        reduction = "umap.harmony", 
        group.by = "seurat_clusters",  # Color the dots by Cluster
        split.by = "condition",        # Split the view by Treatment
        label = TRUE,                  # Show cluster numbers
        ncol = 2) +                    # Arrange side-by-side
  ggtitle("Cluster Split by Condition")
  
  
#Visualize clusters before annotation
#-----------------------------------------------

# Create multi-panel overview
p1 <- DimPlot(seurat_filtered_Humap, reduction = "umap.harmony",
              group.by = "seurat_clusters", label = TRUE, label.size = 5,
              pt.size = 0.1) +
  ggtitle("Clusters (Pre-Annotation)") +
  theme(plot.title = element_text(face = "bold", size = 14))



p2 <- DimPlot(seurat_filtered_Humap, reduction = "umap.harmony",
              group.by = "condition", pt.size = 0.1) +
  ggtitle("Condition") +
  scale_color_manual(values = c("Prerad" = "#2E86AB", 
                                "Rad" = "#F18F01",
                                "Rad6w" = "lightgreen"))
  
##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// find marker genes
###RUN this if using SCTransformer in normalization and clustering
# 1. Switch to standard RNA assay (raw counts)
DefaultAssay(seurat_filtered_Humap) <- "SCT"
# 2. Prep the SCT object for markers (This handles the multiple samples)
# This is the SCT version of 'JoinLayers'
seurat_filtered_Humap <- PrepSCTFindMarkers(seurat_filtered_Humap)
  
# 2. Find markers for every cluster
# min.pct = 0.25: Gene must be in at least 25% of the cluster's cells
# logfc.threshold = 0.25: Gene must be at least slightly upregulated
all_markers <- FindAllMarkers(seurat_filtered_Humap, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
# 3. View the top 10 markers for each cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)  
print(top10_markers, n = 100)

# Save the full list of markers to a CSV file
write.csv(all_markers, "all_cluster_markers_SCT.csv", row.names = FALSE)

# 2. Export this cleaner list
write.csv(top10_markers, "top10_markers_per_cluster_SCT.csv", row.names = FALSE)

saveRDS(seurat_filtered_Humap, file = "preannotated_SCT.rds")
seurat_filtered_Humap <- readRDS("preannotated_SCT.rds")
DimPlot(seurat_filtered_Humap, 
        reduction = "umap.harmony", 
        group.by = "seurat_clusters",  # Color the dots by Cluster
        split.by = "orig.ident",        # Split the view by Treatment
        label = TRUE,                  # Show cluster numbers
        ncol = 2) +                    # Arrange side-by-side
  ggtitle("Cluster Split by Condition")

##
##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// make plot for selective gene in all cluster
# 1. Define the major lineage markers as a named list
major_lineage_markers <- list(
  "T Cells" = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R"),
  "B cells" = c("CD19", "CD38"),
  "Myeloid Cells" = c("LYZ", "CD14", "CD68", "HLA-DRA", "FCGR3A"),
  "Epithelial Cells" = c("EPCAM", "KRT8", "KRT18", "CDH1"),
  "Endothelial Cells" = c("PECAM1", "VWF", "CDH5", "ENG"),
  "Fibroblast" = c("COL1A1", "LUM", "PDPN")
)
# 2. Ensure you are using the RNA assay for visualization
DefaultAssay(seurat_filtered_Humap) <- "RNA"

# 3. Create the DotPlot
# Seurat will automatically group these genes by the names in our list!
DotPlot(seurat_filtered_Humap, features = major_lineage_markers, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  ggtitle("Major Cell Lineages")
##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// merging small clusters /similar clusters together
seurat_filtered_Humap <- RenameIdents(object = seurat_filtered_Humap, 
                                      '7' = '0'
                                     )
print("New Cluster Counts:")
print(table(Idents(seurat_filtered_Humap)))


##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// rename cluster
# 1. Define your new names in order (0, 1, 2...)
# Make sure you have exactly as many names as you have clusters!
new_cluster_ids <- c(
  "Myeloid cells 1",      # Name for Cluster 0
  "CD8",   # Name for Cluster 1
  "CD4",          # Name for Cluster 2
  "Fibroblast",      # Name for Cluster 3
  "Neutrophils",    # Name for Cluster 4
  "B cells",         # Name for Cluster 5
  "Endothelial cells",  # Name for Cluster 6
  "Epithelial cells",
  "Unknown 1", ###cluster 8
  "Myeloid cells 2",
  "Meyloid cells 3",
  "Unkown 2",
  "Unkown 3"
)

# 2. Assign the names to the cluster levels
names(new_cluster_ids) <- levels(seurat_filtered_Humap)

# 3. Rename the identities
seurat_filtered_Humap <- RenameIdents(seurat_filtered_Humap, new_cluster_ids)

# 4. Plot
DimPlot(seurat_filtered_Humap, reduction = "umap.harmony", label = TRUE, pt.size = 0.5) +
  ggtitle("Annotated Cell Types")

DimPlot(seurat_filtered_Humap, 
        reduction = "umap.harmony", 
        split.by = "condition", 
        label = TRUE,
        ncol = 3) + 
  ggtitle("Annotated cell types umap by condition")

# Save the entire object to your hard drive
saveRDS(seurat_filtered_Humap, file = "annotated_seurat_object_SCT.rds")

##load the seurat object
seurat_filtered_Humap <- readRDS("annotated_seurat_object_SCT.rds")

# Verify it worked!
table(seurat_filtered_Humap$condition)
# Verify it looks right
DimPlot(seurat_filtered_Humap, label = TRUE,)

##/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////find out tumor cells
# 1. Create a temporary Seurat object containing ONLY the Prerad condition
seurat_prerad <- subset(seurat_filtered_Humap, subset = condition == "Prerad")

# 2. Plot the UMAP colored by Patient/Sample ID
DimPlot(seurat_prerad, 
        group.by = "orig.ident",   # Ensure this matches your actual Patient/Sample ID column!
        label = FALSE) + 
  ggtitle("UMAP Colored by Patient (Pre-Radiation Only)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 1. Define the classic epithelial and HNSC markers
epithelial_markers <- c("EPCAM", "KRT5", "KRT14")

# 2. Run the FeaturePlot across your entire dataset
FeaturePlot(seurat_filtered_Humap, 
            features = epithelial_markers, 
            ncol = 3,                 # Lines up all 3 plots side-by-side in one row
            cols = c("lightgrey", "red"), # Clean color gradient for visibility
            order = TRUE,             # Pulls the highly expressing cells to the very front
            label = TRUE,             # Keeps the cluster labels on so you can spot the "Unknowns"
            repel = TRUE)

##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// violin plot with p value



# 1. Create your base split Violin Plot
# We save it as an object first to keep the code clean
split_vln <- VlnPlot(seurat_filtered_Humap, 
                     features = "MKI67",           # The gene you want to check
                     split.by = "condition",       # Your Rad vs Prerad column
                     pt.size = 0,
                     cols = c("blue", "red")) +                # Hide the dots for a cleaner look
  ggtitle("FOXP3 Expression: Prerad vs Rad")

# 2. Add the p-values on top!
# The 'aes(group = condition)' tells ggpubr to compare the split halves of each violin
split_vln + stat_compare_means(aes(group = split), 
                               method = "t.test", 
                               label = "p.signif") #

##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// DEG of conditions in cluster

# 1. Create a subset of just Myeloid cells
# Replace "9" with whatever your Myeloid cluster number is
myeloid_cells <- subset(seurat_filtered_Humap, idents = "Myeloid cells")

# 2. (Crucial for v5) Join layers if you haven't already
# This ensures all the data is in one piece for the comparison
myeloid_cells <- JoinLayers(myeloid_cells)

# 1. Set the identity to your treatment column
Idents(myeloid_cells) <- "condition"

# 2. Find genes that change between conditions
# ident.1 = The group of interest (e.g., "Rad")
# ident.2 = The reference group (e.g., "Prerad")
rad_markers <- FindMarkers(myeloid_cells, 
                           ident.1 = "Rad", 
                           ident.2 = "Prerad",
                           min.pct = 0.1,        # Don't need high expression for this
                           logfc.threshold = 0.25) # Look for 25% change or more

# 3. View the top upregulated genes in Radiation
head(rad_markers[order(rad_markers$avg_log2FC, decreasing = TRUE), ])


# Compare a specific gene between Prerad and Rad
VlnPlot(myeloid_cells, features = "CCL22", group.by = "condition", pt.size = 0.1)

# Simple volcano Plot of DEG
# Add a column to identify significant genes (e.g., FC > 1 and p < 0.05)
rad_markers$diffexpressed <- "NO"
rad_markers$diffexpressed[rad_markers$avg_log2FC > 0.5 & rad_markers$p_val_adj < 0.05] <- "UP"
rad_markers$diffexpressed[rad_markers$avg_log2FC < -0.5 & rad_markers$p_val_adj < 0.05] <- "DOWN"

# 2. Create a 'label' column
# We only want to label the top genes, or the plot will be unreadable.
# Let's pick the top 10 genes with the lowest p-value (most significant).
rad_markers$delabel <- NA
top_genes <- head(rad_markers[order(rad_markers$p_val_adj), ], 10)
rad_markers$delabel[row.names(rad_markers) %in% row.names(top_genes)] <- row.names(top_genes)

# Plot it
ggplot(rad_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = delabel)) +
  geom_point(alpha = 0.7) + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  
  # This is the magic part that adds non-overlapping labels
  geom_text_repel(max.overlaps = Inf, box.padding = 0.5) +
  
  labs(title = "Volcano Plot: Radiation vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")


##//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// subcluster myeloid cells

myeloid_cells <- NormalizeData(myeloid_cells)
myeloid_cells <- FindVariableFeatures(myeloid_cells, nfeatures = 2000)
myeloid_cells <- ScaleData(myeloid_cells)

# 2. PCA and UMAP
myeloid_cells <- RunPCA(myeloid_cells)
myeloid_cells <- RunUMAP(myeloid_cells, dims = 1:20)
# 3. Cluster
myeloid_cells <- FindNeighbors(myeloid_cells, dims = 1:20)
myeloid_cells <- FindClusters(myeloid_cells, resolution = 0.1)

# 4. Plot
DimPlot(myeloid_cells, label = TRUE) + ggtitle("Myeloid cells Sub-clusters")
saveRDS(myeloid_cells, file = "Myeloid cells subcluster.rds")

all_markers <- FindAllMarkers(myeloid_cells, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

write.csv(all_markers, "myeloid_cell_subcluster_all_markers.csv", row.names = FALSE)


############////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////GSEA 
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
# 1. Pull the Human Hallmark ("H") gene sets
hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")

# 2. Filter for just IFN Alpha and IFN Gamma
ifn_pathways_df <- hallmark_genes %>%
  filter(gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE"))

# 3. Convert the dataframe into a list format (which fgsea requires)
ifn_pathway_list <- split(x = ifn_pathways_df$gene_symbol, f = ifn_pathways_df$gs_name)

# 1. Run FindMarkers with thresholds set to 0 to get ALL genes
# 1. Manually subset the Seurat object to your target cells first
# Note: Ensure your active Idents() are set to your cell types for this to work
epi_subset <- subset(seurat_filtered_Humap, idents = "CD4")

# 2. Run PrepSCTFindMarkers strictly on the new subset
epi_subset <- PrepSCTFindMarkers(epi_subset)

all_genes_stats <- FindMarkers(epi_subset, 
                               ident.1 = "Rad", 
                               ident.2 = "Prerad", 
                               group.by = "condition", 
                               logfc.threshold = 0,         # Force it to keep everything
                               min.pct = 0,
                               recorrect_umi = FALSE)                 # Force it to keep everything

# 2. Extract the Fold Change and the Gene Names
# fgsea needs a "named numeric vector" sorted from highest to lowest
ranked_genes <- all_genes_stats$avg_log2FC
names(ranked_genes) <- rownames(all_genes_stats)

# 3. Sort the list in decreasing order (Highest logFC at the top)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")

# 2. Filter for just IFN Alpha and IFN Gamma
term2gene_df <- ifn_pathways_df %>% 
  dplyr::select(gs_name, gene_symbol)


# Run clusterProfiler's GSEA
cp_gsea_results <- GSEA(geneList = ranked_genes, 
                        TERM2GENE = term2gene_df, 
                        pvalueCutoff = 0.05,      # Strict cutoff applied here!
                        pAdjustMethod = "BH",     # Benjamini-Hochberg correction
                        minGSSize = 10, 
                        maxGSSize = 500,
                        seed = TRUE)              # Keeps results reproducible


# Plot IFN Alpha
# STEP 1: Pull the exact stats from your results object
# (This looks inside the clusterProfiler object to grab the numbers)
alpha_stats <- cp_gsea_results@result["HALLMARK_INTERFERON_ALPHA_RESPONSE", ]
nes_val <- round(alpha_stats$NES, 2)
padj_val <- signif(alpha_stats$p.adjust, 3)

# Create the custom text (The "\n" puts them on separate lines)
custom_label <- paste0("NES: ", nes_val, "\np.adjust: ", padj_val)


# STEP 2: Generate the plot, but turn OFF the default table
plot_alpha <- gseaplot2(cp_gsea_results, 
                        geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                        title = "IFN Alpha Response in CD4 (Rad vs Prerad)", 
                        pvalue_table = FALSE,   # <--- Turned off!
                        base_size = 14)


# STEP 3: Inject your custom text directly into the top plot
# (plot_alpha[[1]] specifically targets the green curve section at the top)
plot_alpha[[1]] <- plot_alpha[[1]] + 
  annotate("text", 
           x = Inf, y = Inf,        # Automatically finds the top-right corner
           label = custom_label, 
           hjust = 1.2, vjust = 2,  # Nudges the text slightly inward so it isn't cut off
           size = 8, 
           fontface = "bold",
           color = "black")

# View the beautifully customized plot!
plot_alpha
# Plot IFN Gamma
# STEP 1: Pull the exact stats for IFN Gamma
gamma_stats <- cp_gsea_results@result["HALLMARK_INTERFERON_GAMMA_RESPONSE", ]
nes_val_gamma <- round(gamma_stats$NES, 2)
padj_val_gamma <- signif(gamma_stats$p.adjust, 3)

# Create the custom text
custom_label_gamma <- paste0("NES: ", nes_val_gamma, "\np.adj: ", padj_val_gamma)


# STEP 2: Generate the base plot, with the default table turned OFF
plot_gamma <- gseaplot2(cp_gsea_results, 
                        geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                        title = "IFN Gamma Response in CD4 (Rad vs Prerad)", 
                        pvalue_table = FALSE,   
                        base_size = 14)


# STEP 3: Inject your custom text into the top-right corner of the curve
plot_gamma[[1]] <- plot_gamma[[1]] + 
  annotate("text", 
           x = Inf, y = Inf,        
           label = custom_label_gamma, 
           hjust = 1.2, vjust = 2,  
           size = 8, 
           fontface = "bold",
           color = "black")

# View the beautifully customized IFN Gamma plot!
plot_gamma


###/////////////////////////////////////////////////////////////////////////////////////////////////////////////////double check the ifn socre with addmodulescore
# 1. Pull out the IFN-Gamma genes into a standard vector
ifng_genes <- ifn_pathway_list[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]
ifna_genes <- ifn_pathway_list[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]

# 2. Put that vector inside a list() - THIS IS REQUIRED FOR SEURAT!
genes_to_score <- list(ifng_genes)

# Run the scoring algorithm
seurat_filtered_Humap <- AddModuleScore(
  object = seurat_filtered_Humap,
  features = genes_to_score,
  name = "IFNg_Score"  # What you want to call the new column
)
summary(seurat_filtered_Humap$IFNg_Score1)
# Plot the module score across all clusters, split by condition
VlnPlot(seurat_filtered_Humap, 
        features = "IFNg_Score1",        # Use the new score column
        group.by = "seurat_clusters",    # Group by your clusters on the X-axis
        split.by = "condition",          # Split the violins by Rad vs Prerad
        pt.size = 0) +                   # Hide the dots for a clean look
  
  # Add a clean title
  ggtitle("IFN-GAMMA Module Score by Cluster") +
  
  # Optional: Add significance stars just like you did earlier!
  stat_compare_means(aes(group = split), 
                     method = "wilcox.test", 
                     label = "p.signif") +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
library(patchwork)
# Note: We removed the 'cols' argument from here because ggplot will handle it below
# Generate the plot without the 'cols' argument
FeaturePlot(seurat_filtered_Humap, 
            features = "IFNg_Score1", 
            split.by = "condition", 
            order = TRUE,
            label = TRUE,           
            repel = TRUE) +
  
  # Add the Master Title
  plot_annotation(title = "IFN-Gamma Pathway Enrichment)") &
  
  # Create the zero-centered diverging color bar
  scale_color_gradient2(
    low = "steelblue",       # Negative values (suppressed)
    mid = "lightgreen",           # Zero value (background level)
    high = "red",     # Positive values (activated)
    midpoint = 0,            # Locks the white color exactly at a score of 0
    name = "Value"           # Or "Z-Score" / "Module Score"
  ) &
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
# Define your custom baseline cutoff based on your summary stats
# (e.g., if your median score is 0.05, set this to 0.05)
my_cutoff <- 0.05 
max_score <- max(seurat_filtered_Humap$IFNg_Score1)

FeaturePlot(seurat_filtered_Humap, 
            features = "IFNg_Score1", 
            split.by = "condition", 
            order = TRUE,
            label = TRUE,           
            repel = TRUE) +
  
  plot_annotation(title = "IFN-Gamma Pathway Enrichment") &
  
  # Use a 2-color gradient to make contrast pop
  scale_color_gradient(
    low = "lightgrey",             # Everything at or below your cutoff becomes lightgrey
    high = "blue",                  # High expression stays vibrant red
    limits = c(my_cutoff, max_score), 
    oob = scales::squish,          # <--- Squishes all cells below 'my_cutoff' to the 'low' color
    name = "Score"
  ) &
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))


#####//////////////////////////////////////////////////////////////////////////////GSEA for all clusters in dot plot
# 1. Get a list of all the cluster names currently in your object
all_clusters <- levels(Idents(seurat_filtered_Humap))

# 2. Create an empty list to store the results
gsea_summary_list <- list()

# 3. Loop through every single cluster
for (clust in all_clusters) {
  
  # try() prevents the loop from crashing if a rare cluster lacks Rad/Prerad cells
  try({
    # A. Run FindMarkers for ALL genes in this specific cluster
    stats <- FindMarkers(seurat_filtered_Humap, 
                         ident.1 = "Rad", 
                         ident.2 = "Prerad", 
                         group.by = "condition", 
                         subset.ident = clust,
                         logfc.threshold = 0, 
                         min.pct = 0,
                         verbose = FALSE)
    
    # B. Rank the genes by log2 Fold Change
    ranked_genes <- stats$avg_log2FC
    names(ranked_genes) <- rownames(stats)
    ranked_genes <- sort(ranked_genes, decreasing = TRUE)
    
    # C. Run fgsea against our IFN pathway list
    res <- fgsea(pathways = ifn_pathway_list, 
                 stats = ranked_genes,
                 minSize = 10, maxSize = 500)
    
    # D. Clean up the results and add the Cluster name
    # (We drop the "leadingEdge" column because it breaks the merge step later)
    res_clean <- res %>% 
      dplyr::select(pathway, pval, padj, NES) %>% 
      mutate(Cluster = clust)
    
    # Save into our list
    gsea_summary_list[[clust]] <- res_clean
    
  }, silent = TRUE) # End of try block
}

# 4. Stitch all the individual cluster results into one master dataframe
final_gsea_df <- bind_rows(gsea_summary_list)


# Transform the Adjusted P-value so smaller numbers make bigger dots (-log10)
final_gsea_df <- final_gsea_df %>%
  mutate(Significance = -log10(padj))

# Clean up the pathway names so they look nice on the plot
final_gsea_df$pathway <- gsub("HALLMARK_INTERFERON_", "IFN-", final_gsea_df$pathway)
final_gsea_df$pathway <- gsub("_RESPONSE", "", final_gsea_df$pathway)

# Generate the Summary Dot Plot
ggplot(final_gsea_df, aes(x = Cluster, y = pathway, color = NES, size = Significance)) +
  geom_point() +
  
  # Create a beautiful Red/Blue color scale centered at 0
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                        name = "NES\n(Red = Up in Rad)") +
  
  # Set the dot sizes
  scale_size_continuous(range = c(2, 10), name = "-log10(p.adj)") +
  
  # Clean aesthetics
  theme_bw() +
  labs(title = "Interferon Pathway Enrichment Across All Clusters",
       subtitle = "Comparing Rad vs Prerad",
       x = "Cell Cluster",
       y = "Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12))


