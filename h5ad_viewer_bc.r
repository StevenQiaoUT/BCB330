if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)

p <- path.expand("~/Desktop/UofT Files/Third Year/BCB330/at.h5ad")
p <- normalizePath(p, mustWork = TRUE)

sce <- readH5AD(p)

sce
colnames(colData(sce))   
rownames(rowData(sce))   
assayNames(sce)          

# view metadata about column so cells
head(colData(sce))

# summarize counts of each categorical variable within CellType key
table(colData(sce)$CellType)

# view metadata about rows so genes
head(rowData(sce))

# summarize counts of the booleans if gene is highly variable
table(rowData(sce)$highly_variable)

# express matrix and use "X" to determine matrix type
exprs <- assay(sce, "X")

dim(exprs)

exprs[1:5, 1:5]

exprs[1:20, 1:20]

# See the clustering results from leiden

# Shows metadata/stored results and analysis of the entire .h5ad file 
names(metadata(sce))  

metadata(sce)$leiden  

# Conducted PCA reduction to two dimensions
reducedDimNames(sce)  
umap_coords <- reducedDim(sce, "X_umap_hm")
head(umap_coords)       

plot_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  color = colData(sce)[["CellType"]]
)

# plot by cell type
library(ggplot2)
ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = color)) +
  geom_point(size = 0.6, alpha = 0.8) +
  coord_equal() +
  labs(title = "UMAP of cells", color = "CellType") +
  theme_minimal()

# plot by gene intensity, AP006222.2 chosen for now
library(scales)

# Choose a gene
gene <- "AP006222.2"

# See all the different types of matrices
assayNames(sce)

# Choose the first matrix
assay_name <- assayNames(sce)[1]

# Store the expression levels into a vector 
# Pull up the assay first
# In our file, we only have 1 matrix "X"
# In this code, [gene, ] syntax is ["AP006222.2", ], so we are taking all the 
# expression levels from the AP006222.2 column
expr_vec <- as.numeric(assay(sce, assay_name)[gene, ])

# Data points and properties we want to plot
plot_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  expr  = expr_vec
)

# plotting nitty-gritties
ggplot(plot_df, aes(UMAP1, UMAP2, color = expr)) +
  geom_point(size = 0.6, alpha = 0.9) +
  coord_equal() +
  scale_color_gradient(
    low = "grey70",     # 0 expression → gray
    high = "red",       # high expression → deep red
    name = gene
  ) +
  labs(title = paste("UMAP colored by", gene, "expression")) +
  theme_minimal()


# Packages
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("zellkonverter", quietly = TRUE)) BiocManager::install("zellkonverter")
if (!requireNamespace("bench", quietly = TRUE)) install.packages("bench")

library(reticulate)
library(zellkonverter)
library(bench)

# Paths (edit these)
zarr_path <- normalizePath("~/Desktop/UofT Files/Third Year/BCB330/new_at.zarr", mustWork = TRUE)  # folder
h5ad_path <- normalizePath("~/Desktop/UofT Files/Third Year/BCB330/at.h5ad",       mustWork = TRUE)

# Import anndata once so we don't include import overhead in timings
anndata <- import("anndata", convert = FALSE)


load_zarr_as_sce <- function(path) {
  ad <- anndata$read_zarr(path)            # Python read
  zellkonverter::AnnData2SCE(ad)           # convert to SCE
}

load_h5ad_as_sce <- function(path) {
  zellkonverter::readH5AD(path)            # direct to SCE
}

# Warm-up (optional) to avoid first-call overhead in the benchmark
invisible(try(load_zarr_as_sce(zarr_path), silent = TRUE))
invisible(try(load_h5ad_as_sce(h5ad_path), silent = TRUE))
gc()

res_sce <- bench::mark(
  zarr_to_sce = load_zarr_as_sce(zarr_path),
  h5ad_to_sce = load_h5ad_as_sce(h5ad_path),
  iterations = 5,
  check = FALSE,       # we won’t compare object equality here
  filter_gc = TRUE
)

print(res_sce[, c("expression","min","median","itr/sec","mem_alloc","gc")])


