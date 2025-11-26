library(reticulate)
library(zellkonverter)
library(bench)

zarr_path <- normalizePath("~/Desktop/UofT Files/Third Year/BCB330/new_at.zarr", mustWork = TRUE)
h5ad_path <- normalizePath("~/Desktop/UofT Files/Third Year/BCB330/at.h5ad",    mustWork = TRUE)

anndata <- import("anndata", convert = FALSE)

# Functions that return NULL so bench doesn't store giant results
zarr_to_sce_null <- function() {
  ad  <- anndata$read_zarr(zarr_path)
  sce <- zellkonverter::AnnData2SCE(ad)
  rm(ad, sce); gc()
  NULL
}
h5ad_to_sce_null <- function() {
  sce <- zellkonverter::readH5AD(h5ad_path)
  rm(sce); gc()
  NULL
}

# Warm-up to avoid first-call overhead
invisible(zarr_to_sce_null()); invisible(h5ad_to_sce_null()); gc()

res_sce <- bench::mark(
  zarr_to_sce = zarr_to_sce_null(),
  h5ad_to_sce = h5ad_to_sce_null(),
  iterations  = 3,           # increase if you have headroom
  check       = FALSE,
  filter_gc   = TRUE
)

print(res_sce[, c("expression","min","median","itr/sec","mem_alloc","gc")])
