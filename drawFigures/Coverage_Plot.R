library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)

# We need the reticulate package to interact with Python and read h5ad files
library(reticulate)
# Import anndata through reticulate
anndata <- import("anndata")

counts <- Read10X_h5("D:/Desktop/R/SignacMotif/data/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "D:/Desktop/R/SignacMotif/data/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = 'D:/Desktop/R/SignacMotif/data/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

# Read AnnData file
adata <- anndata$read_h5ad("D:/Desktop/R/SignacMotif/orig.h5ad")

# Extract cell type labels
cell_types <- py_to_r(adata$obs["Pred"])
head(cell_types)

# Ensure that row names in cell_types match cell barcodes in Seurat object
# You may need to adjust barcode format to match between the two objects
head(colnames(brain))
head(rownames(cell_types))

# Add cell type labels to Seurat object
brain$cell_type <- cell_types[match(colnames(brain), rownames(cell_types)), "Pred"]

# Verify that labels were added correctly
table(brain$cell_type, useNA = "ifany")

# Add gene annotations to Seurat object
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'  # Convert chromosome naming to UCSC style (chr1, chr2, etc.)

# Check consistency of annotation chromosome naming style with peaks assay
genome(annotations) <- "mm10"
Annotation(brain) <- annotations

# Now try using gene names
DefaultAssay(brain) <- "peaks"

# Neurod6  Tmem119  Ndrg2  Mag  Itm2a  Gad2  Sst
# First, make sure you have the required font
library(svglite)
library(extrafont)
# You may need to import fonts first if you haven't done so
# font_import()
# loadfonts()

p <- CoveragePlot(
  object = brain,
  region = "Gad2",
  extend.upstream = 3000,
  extend.downstream = 3000,
  group.by = "cell_type",
  heights = c(8, 1),
  peaks = FALSE,
  annotation = "gene"
)

# Custom theme settings
custom_theme <- theme(
  # Basic text elements
  text = element_text(family = "Times New Roman", size = 12),
  axis.text = element_text(family = "Times New Roman", size = 8),
  axis.title = element_text(family = "Times New Roman", size = 12),
  legend.text = element_text(family = "Times New Roman", size = 14),
  legend.title = element_text(family = "Times New Roman", size = 14),
  axis.title.x = element_blank()
)

# Apply theme
p <- p & custom_theme

# # Save the plot as SVG
# ggsave(
#   filename = "D:/Desktop/R/SignacMotif/Coverage_plot/Gad2_coverage_plot.svg",
#   plot = p,
#   device = "svg",
#   width = 4,
#   height = 5
# )