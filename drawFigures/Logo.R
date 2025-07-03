# First, check if you have saved the necessary objects
# If you've already run Signac.R, you can load the saved objects:
seurat_obj <- readRDS("D:/Desktop/R/SignacMotif/seurat_signac_object.rds")
motif_results <- readRDS("D:/Desktop/R/SignacMotif/motif_enrichment_results.rds")

# Load required libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(svglite)
library(extrafont)

# Load fonts
loadfonts()

# Now run the plot_top_motifs function from section 5.1
plot_top_motifs <- function() {
  # Ensure svg save functionality is available
  if(!require(svglite)) {
    install.packages("svglite")
    library(svglite)
  }
  
  # For each cell type
  for(cell_type in names(motif_results)) {
    # Get the top 3 most significant motifs
    top_motifs <- head(rownames(motif_results[[cell_type]][order(motif_results[[cell_type]]$pvalue), ]), 3)
    
    if(length(top_motifs) > 0) {
      # Create motif logos using default font
      p <- MotifPlot(
        object = seurat_obj,
        motifs = top_motifs,
        assay = "ATAC"
      )
      
      # Modify plot appearance: set x-axis labels vertical, enlarge and bold motif names
      p <- p + theme(
        text = element_text(family = "Times New Roman"),  # Set font to Times New Roman
        axis.text.x = element_text(size = 17, angle = 90, hjust = 1), # Make x-axis labels vertical and bold
        axis.text.y = element_text(size = 18), # Set y-axis label size and bold
        axis.title.x = element_text(size = 21),  # Set x-axis title font size and bold
        axis.title.y = element_text(size = 24, margin = margin(r = 10)),  # Set y-axis title font size and bold
        plot.title = element_text(size = 15),   # Set plot title font size and bold
        axis.ticks.length = unit(0.2, "cm"),  # Set axis tick mark length
        strip.text = element_text(size = 24)  # Enlarge motif names and make bold
      )
      
      # Save as SVG
      svg_filename <- paste0("D:/Desktop/R/SignacMotif/", cell_type, "_motif_logos.svg")
      svglite(svg_filename, width = 11, height = 4)
      print(p)
      dev.off()
    }
  }
  
  return(TRUE)
}

# Execute the function
cat("Generating motif plots with Signac...\n")
plot_result <- plot_top_motifs()
cat("Completed generating motif plots\n")