# Single-cell ATAC-seq Motif Analysis Pipeline
# This script performs motif enrichment analysis on cell-type specific peaks

# ========================================
# 1. LOAD REQUIRED PACKAGES
# ========================================

# Core packages for single-cell ATAC-seq analysis
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)            # mm10 annotation
library(BSgenome.Mmusculus.UCSC.mm10)   # mm10 reference genome

# Data manipulation and visualization
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(hdf5r)
library(anndata)
library(reticulate)

# Parallel processing
library(parallel)
library(doParallel)
library(foreach)

# Motif analysis
library(TFBSTools)
library(JASPAR2020)

# Font and graphics support
library(extrafont)
library(svglite)

# ========================================
# 2. SETUP PARALLEL PROCESSING
# ========================================

# Import system fonts (run once for first time)
# font_import()
loadfonts()
theme_set(theme_minimal()) 

# Set up parallel processing (leave one core for system)
cores <- detectCores() - 1
registerDoParallel(cores)

# Set memory limit for parallel processing (24GB)
options(future.globals.maxSize = 24 * 1024^3)

# ========================================
# 3. DATA IMPORT AND PREPROCESSING
# ========================================

# Load data from h5ad file
# Convert Python AnnData object to R Seurat object
ad <- read_h5ad("D:/Desktop/R/SignacMotif/orig.h5ad")

# Extract necessary information
counts <- t(as(ad$X, "dgCMatrix"))       # Transpose sparse matrix
cell_metadata <- as.data.frame(ad$obs)   # Cell metadata
peak_metadata <- as.data.frame(ad$var)   # Peak metadata

# Create initial Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  assay = "ATAC",
  meta.data = cell_metadata
)

# ========================================
# 4. BUILD SIGNAC CHROMATIN ASSAY
# ========================================

# Create GRanges object from peak metadata
peaks_gr <- GRanges(
  seqnames = peak_metadata$chr,
  ranges = IRanges(start = peak_metadata$start, end = peak_metadata$end)
)

# Get chromosome length information
seqinfo <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
genome_name <- "mm10"

# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = genome_name,
  ranges = peaks_gr,
  fragments = NULL  # Add fragments file if available
)

# Replace assay in Seurat object
seurat_obj[["ATAC"]] <- chrom_assay

# ========================================
# 5. FIND CELL-TYPE SPECIFIC PEAKS (1K VERSION)
# ========================================

# Set cell identity using Pred column
Idents(seurat_obj) <- "Pred"
cell_types <- unique(seurat_obj$Pred)

# Storage for differential accessibility peaks
diff_peaks_list <- list()

# Find top 1000 specific peaks for each cell type
for (cell_type in cell_types) {
  cat("Processing cell type:", cell_type, "\n")
  
  # Find differential peaks using likelihood ratio test
  da_peaks <- FindMarkers(
    object = seurat_obj,
    ident.1 = cell_type,
    ident.2 = NULL,
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'n_genes'
  )
  
  # Sort by p-value
  da_peaks <- da_peaks[order(da_peaks$p_val), ]
  
  # Select top 1000 specific peaks (or all available if < 1000)
  if (nrow(da_peaks) >= 1000) {
    top_da_peaks <- rownames(head(da_peaks, 1000))
  } else {
    top_da_peaks <- rownames(da_peaks)
    cat("Warning: Only", length(top_da_peaks), "peaks available for", cell_type, "\n")
  }
  
  diff_peaks_list[[cell_type]] <- top_da_peaks
}

# Add names to results list
names(diff_peaks_list) <- cell_types

# Print peak counts for each cell type
peak_counts <- sapply(diff_peaks_list, length)
print(peak_counts)

# Save cell-type specific peaks list (1K version)
saveRDS(diff_peaks_list, "D:/Desktop/R/SignacMotif/cell_type_specific_peaks.rds")

# ========================================
# 6. FIND CELL-TYPE SPECIFIC PEAKS (10K VERSION)
# ========================================

# Storage for 10K version
diff_peaks_list_10k <- list()

# Find top 10000 specific peaks for each cell type
for (cell_type in cell_types) {
  cat("Processing cell type:", cell_type, " (10k version)\n")
  
  # Find differential peaks
  da_peaks <- FindMarkers(
    object = seurat_obj,
    ident.1 = cell_type,
    ident.2 = NULL,
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'n_genes'
  )
  
  # Sort by p-value
  da_peaks <- da_peaks[order(da_peaks$p_val), ]
  
  # Select top 10000 specific peaks (or all available if < 10000)
  if (nrow(da_peaks) >= 10000) {
    top_da_peaks <- rownames(head(da_peaks, 10000))
  } else {
    top_da_peaks <- rownames(da_peaks)
    cat("Warning: Only", length(top_da_peaks), "peaks available for", cell_type, "\n")
  }
  
  diff_peaks_list_10k[[cell_type]] <- top_da_peaks
}

# Add names to results list
names(diff_peaks_list_10k) <- cell_types

# Print peak counts for each cell type (10K version)
peak_counts_10k <- sapply(diff_peaks_list_10k, length)
print(peak_counts_10k)

# Save cell-type specific peaks list (10K version)
saveRDS(diff_peaks_list_10k, "D:/Desktop/R/SignacMotif/cell_type_specific_peaks_10k.rds")

# ========================================
# 7. MOTIF ENRICHMENT ANALYSIS
# ========================================

# Get JASPAR PWMs (Position Weight Matrices)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information to Seurat object
seurat_obj <- AddMotifs(
  object = seurat_obj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# Storage for motif enrichment results
motif_results <- list()

# Perform motif enrichment analysis for each cell type
for(cell_type in names(diff_peaks_list)) {
  cat("Processing motif enrichment for cell type:", cell_type, "\n")
  features <- diff_peaks_list[[cell_type]]
  
  # Skip if not enough differential peaks
  if(length(features) < 10) {
    motif_results[[cell_type]] <- NULL
    next
  }
  
  # Perform motif enrichment analysis
  motif_results[[cell_type]] <- FindMotifs(
    object = seurat_obj,
    features = features,
    pwm = pfm,
    background.size = 40000
  )
}

# Remove NULL results
motif_results <- motif_results[!sapply(motif_results, is.null)]

# ========================================
# 8. VISUALIZATION FUNCTIONS
# ========================================

# Function to create motif logo plots
plot_top_motifs <- function() {
  # Ensure svg support is available
  if(!require(svglite)) {
    install.packages("svglite")
    library(svglite)
  }
  
  # For each cell type
  for(cell_type in names(motif_results)) {
    # Get top 3 most significant motifs
    top_motifs <- head(rownames(motif_results[[cell_type]][order(motif_results[[cell_type]]$pvalue), ]), 3)
    
    if(length(top_motifs) > 0) {
      # Create motif logo plot
      p <- MotifPlot(
        object = seurat_obj,
        motifs = top_motifs,
        assay = "ATAC"
      )
      
      # Customize plot appearance
      p <- p + theme(
        text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"),
        strip.text = element_text(size = 20)
      )
      
      # Save as SVG (COMMENTED OUT)
      # svg_filename <- paste0("D:/Desktop/R/SignacMotif/", cell_type, "_motif_logos.svg")
      # svglite(svg_filename, width = 8, height = 4)
      # print(p)
      # dev.off()
    }
  }
  
  return(TRUE)
}

# Function to prepare UpSet plot data
prepare_upset_data <- function() {
  # Select top 40 significant motifs for each cell type
  top_motifs_by_cell <- list()
  
  for(cell_type in names(motif_results)) {
    results <- motif_results[[cell_type]]
    # Sort by adjusted p-value
    results <- results[order(results$p.adjust), ]
    # Select top 40
    top_motifs <- head(rownames(results), 40)
    top_motifs_by_cell[[cell_type]] <- top_motifs
  }
  
  # Create shared motif matrix
  all_motifs <- unique(unlist(top_motifs_by_cell))
  motif_matrix <- matrix(0, nrow = length(all_motifs), ncol = length(top_motifs_by_cell))
  rownames(motif_matrix) <- all_motifs
  colnames(motif_matrix) <- names(top_motifs_by_cell)
  
  for(i in 1:length(top_motifs_by_cell)) {
    cell_type <- names(top_motifs_by_cell)[i]
    motifs <- top_motifs_by_cell[[cell_type]]
    motif_matrix[motifs, cell_type] <- 1
  }
  
  return(list(
    top_motifs_by_cell = top_motifs_by_cell,
    motif_matrix = motif_matrix
  ))
}

# ========================================
# 9. GENERATE PLOTS (COMMENTED OUT)
# ========================================

# Generate motif plots
cat("Generating motif plots with Signac...\n")
plot_result <- plot_top_motifs()

# Prepare UpSet plot data
upset_data <- prepare_upset_data()

# Create and save UpSet plot (COMMENTED OUT)
# if(require(UpSetR)) {
#   # Get available set names
#   available_sets <- colnames(upset_data$motif_matrix)
#   
#   # Remove specific cell types if needed
#   sets_to_remove <- c("Collisions", "Astrocytes")
#   sets_to_use <- setdiff(available_sets, sets_to_remove)
#   
#   # Check if enough sets remain
#   if(length(sets_to_use) < 2) {
#     warning("Not enough sets remaining after removal for meaningful UpSet plot")
#     sets_to_use <- available_sets
#   }
#   
#   # Reverse set order for diagonal appearance on left
#   reversed_sets <- rev(sets_to_use)
#   
#   # Create filtered data frame
#   upset_data_filtered <- as.data.frame(upset_data$motif_matrix[, sets_to_use])
#   
#   # Create UpSet plot
#   upset_plot <- upset(
#     upset_data_filtered,
#     sets = reversed_sets,
#     mb.ratio = c(0.6, 0.4),
#     order.by = "degree",
#     decreasing = FALSE,
#     group.by = "degree",
#     keep.order = TRUE,
#     mainbar.y.label = "Intersection Size",
#     text.scale = 1.2,
#     sets.x.label = "Set Size"
#   )
#   
#   # Save plot (COMMENTED OUT)
#   # pdf("D:/Desktop/R/SignacMotif/upset_plot.pdf", width = 12, height = 8)
#   # print(upset_plot)
#   # dev.off()
# }

# ========================================
# 10. SAVE RESULTS
# ========================================

# Save processed Seurat object
saveRDS(seurat_obj, "D:/Desktop/R/SignacMotif/seurat_signac_object.rds")

# Save motif enrichment results
saveRDS(motif_results, "D:/Desktop/R/SignacMotif/motif_enrichment_results.rds")

# Clean up parallel processing
stopImplicitCluster()

cat("Analysis completed successfully!\n")