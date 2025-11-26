# ===== 0) Clean + packages =====
# Restart R first (Session -> Restart R) before sourcing this file.

# Remove env vars that can force the wrong Python
Sys.unsetenv(c("RETICULATE_PYTHON", "VIRTUAL_ENV", "CONDA_PREFIX", "PYTHONPATH"))

if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("zellkonverter", quietly = TRUE)) BiocManager::install("zellkonverter")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(reticulate)
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)

# ===== 1) Ensure Miniconda and pin the conda env BEFORE Python initializes =====
if (is.null(reticulate::conda_binary())) {
  reticulate::install_miniconda()
}

# Create env once if needed
env_name <- "r-anndata"
envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
if (!(env_name %in% envs)) {
  reticulate::conda_create(env_name, packages = c("python=3.11"))
}

# Hard-pin Python by path and tell reticulate before initialization:
conda_py <- reticulate::conda_python(env_name)   # e.g., ~/Library/r-miniconda-*/envs/r-anndata/bin/python
Sys.setenv(RETICULATE_PYTHON = conda_py)

# Do NOT call py_config() yet — that would initialize Python too early.

# Install Python deps into this env (safe to re-run)
reticulate::conda_install(env_name, c("numpy", "scipy"))
reticulate::py_install(c("anndata>=0.9", "zarr"), envname = env_name, pip = TRUE)

# Now initialize and confirm the chosen Python
reticulate::py_config()

# ===== 2) Read .zarr and convert to SCE =====
anndata <- import("anndata", convert = FALSE)

zarr_path <- normalizePath(path.expand("~/Desktop/UofT Files/Third Year/BCB330/new_at.zarr"), mustWork = TRUE)
adata <- anndata$read_zarr(zarr_path)
sce   <- zellkonverter::AnnData2SCE(adata)

# Optional peek
sce
SummarizedExperiment::assayNames(sce)
SingleCellExperiment::reducedDimNames(sce)

# ===== 3) Plot UMAP colored by a gene =====
# Your reductions are: "X_pca", "X_tsne", "X_umap"
umap_coords <- SingleCellExperiment::reducedDim(sce, "X_umap")
stopifnot(nrow(umap_coords) == ncol(sce))

# Expression matrix: choose an assay; you have "X" and "RNA"
assay_name <- "X"  # change to "RNA" if that’s what you want
mat <- SummarizedExperiment::assay(sce, assay_name)

gene <- "ATTLL1"  # change as needed
g <- match(gene, rownames(sce))
if (is.na(g)) {
  candidates <- grep(gene, rownames(sce), value = TRUE, ignore.case = TRUE)
  stop(sprintf("Gene '%s' not found.\nTop matches: %s", gene, paste(head(candidates, 10), collapse = ", ")))
}

# Sparse check: use class-based test (Matrix::isSparse is not exported)
is_sparse <- inherits(mat, "sparseMatrix")
expr_vec  <- if (is_sparse) as.numeric(mat[g, ]) else as.numeric(mat[g, ])

plot_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  expr  = expr_vec
)

ggplot(plot_df, aes(UMAP1, UMAP2, color = expr)) +
  geom_point(size = 0.5, alpha = 0.9) +
  coord_equal() +
  scale_color_gradient(low = "grey70", high = "red", name = gene) +
  labs(title = paste("UMAP colored by", gene, "expression (assay:", assay_name, ")")) +
  theme_minimal()

