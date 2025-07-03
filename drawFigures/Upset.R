# Load necessary packages
library(UpSetR)  # For UpSet plots
library(parallel)  # For parallel processing
library(extrafont)  # For font handling
library(svglite)  

# Load fonts
loadfonts()

# Try to set Windows fonts
tryCatch({
  windowsFonts(
    times = windowsFont("Times New Roman")
  )
  cat("Successfully loaded Windows fonts\n")
}, error = function(e) {
  cat("Failed to load Windows fonts:", e$message, "\n")
})

par(family = "times")

# Load previously saved results
seurat_obj <- readRDS("D:/Desktop/R/SignacMotif/seurat_signac_object.rds")
motif_results <- readRDS("D:/Desktop/R/SignacMotif/motif_enrichment_results.rds")

# 5.2 Prepare data for UpSet plot display
prepare_upset_data <- function() {
  # Select top 40 significant motifs for each cell type
  top_motifs_by_cell <- list()
  
  for(cell_type in names(motif_results)) {
    results <- motif_results[[cell_type]]
    # Sort by p-value
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

# Run UpSet analysis
upset_data <- prepare_upset_data()

# Create and save UpSet plot
# Get set names from data frame
available_sets <- colnames(upset_data$motif_matrix)

# Remove "Collisions" and "Astrocytes" cell types (if they exist)
sets_to_remove <- c("Collisions", "Astrocytes")
sets_to_use <- setdiff(available_sets, sets_to_remove)

# If not enough sets remain, issue a warning
if(length(sets_to_use) < 2) {
  warning("Insufficient sets remaining after removing cell types, cannot create meaningful UpSet plot")
  # Still continue, using all available sets
  sets_to_use <- available_sets
}

# Reverse set order to make diagonal appear on the left
reversed_sets <- rev(sets_to_use)

# Create data frame using only selected sets
upset_data_filtered <- as.data.frame(upset_data$motif_matrix[, sets_to_use])

# Create UpSet plot
upset_plot <- upset(
  upset_data_filtered,
  sets = reversed_sets,
  mb.ratio = c(0.6, 0.4),
  order.by = "degree",
  decreasing = FALSE,
  group.by = "degree",
  keep.order = TRUE,
  mainbar.y.label = "Intersection Size",
  text.scale = 2.5,
  sets.x.label = "set size",
  point.size = 3.0,  # Increase point size
  line.size = 1.1    # Increase line thickness
)

# Save as PNG format
png("D:/Desktop/R/SignacMotif/upset_plot.png", width = 12*1000, height = 8*1000, res = 1000, family = "times")
print(upset_plot)
dev.off()
cat("PNG saved using times font\n")

# Save SVG using svglite
svglite("D:/Desktop/R/SignacMotif/upset_plot.svg", width = 12, height = 8,
        system_fonts = list(sans = "Times New Roman", serif = "Times New Roman", mono = "Times New Roman"))
print(upset_plot)
dev.off()

# Close parallel cluster and release resources
stopImplicitCluster()