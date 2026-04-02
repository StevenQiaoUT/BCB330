# Get unique cell types
unique_celltypes <- unique(cell_types)

print(cell_types)
print(unique_celltypes)

# Initialize matrix to store results: genes (rows) x cell types (columns)
avg_expr_matrix <- matrix(
  nrow = nrow(exprs),
  ncol = length(unique_celltypes),
  dimnames = list(rownames(exprs), unique_celltypes)
)

# For each cell type, compute mean expression
for (ct in unique_celltypes) {
  # Find cells of this type
  cell_indices <- which(cell_types == ct)
  
  # Calculate mean expression across those cells for each gene
  if (length(cell_indices) > 1) {
    avg_expr_matrix[, ct] <- rowMeans(exprs[, cell_indices], na.rm = TRUE)
  } else {
    avg_expr_matrix[, ct] <- exprs[, cell_indices]
  }
}

library(pheatmap)
# View result
head(avg_expr_matrix)
dim(avg_expr_matrix)  # genes x cell_types

expressed <- rowSums(avg_expr_matrix) > 0.01
pheatmap(log1p(avg_expr_matrix[expressed, ][1:50, ]),
         main = "Average expression by cell type (top 50 genes)",
         fontsize_row = 6)
