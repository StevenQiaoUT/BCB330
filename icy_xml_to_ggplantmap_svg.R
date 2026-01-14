library(ggplot2)
library(ggPlantmap)
library(tidyverse)

setwd("/Users/stevenqiao/Desktop/ggplantmap_data")

# Load data
data <- read.csv("avg_expression.csv")
reference <- read.table("celltype_reference.txt", sep="\t", header=TRUE)

print(paste("Loaded", nrow(data), "genes"))
print(paste("Loaded", nrow(reference), "cell types"))

# Create mapping from your cell types to ggPlantmap ROI names
celltype_mapping <- tribble(
  ~Identity, ~ROI_mapped,
  "0: Mesophyll_1", "Parenchima.palisade",
  "1: Mesophyll_2", "Parenchima.palisade",
  "2: Mesophyll_3", "Parenchima.palisade",
  "3: Mesophyll_4", "Parenchima.palisade",
  "4: Mesophyll_5", "Parenchima.sponge",
  "5: Epidermal_1", "epidermis.adaxial",
  "6: Immune active", "Parenchima.sponge",
  "7: Mesophyll_6", "Parenchima.sponge",
  "8: Guard_1", "epidermis.stomata",
  "9: Mesophyll_7", "Parenchima.sponge",
  "10: Defense state", "Parenchima.sponge",
  "11: Vascular_1", "vascularbundle.xylem",
  "12: Vascular_2", "vascularbundle.phloem",
  "13: Epidermal_2", "epidermis.abaxial",
  "14: Phloem companion", "vascularbundle.phloem",
  "15: Stress responsive", "Parenchima.sponge",
  "16: Phloem Parenchyma", "vascularbundle.phloem",
  "17: Sieve element_responsive", "vascularbundle.phloem",
  "18: Dividing_1", "Parenchima.palisade",
  "19: Guard_2", "epidermis.stomata",
  "20: Dividing_2", "Parenchima.sponge",
  "21: Epidermal_3", "epidermis.abaxial",
  "22: Metabolic stress state", "Parenchima.sponge",
  "23: Hydathode", "epidermis.abaxial",
  "24: Trichome", "epidermis.adaxial",
  "25: Myrosin", "Parenchima.palisade",
  "26: Sugar metabolic state", "Parenchima.sponge"
)

# Merge reference with mapping
reference <- reference %>%
  left_join(celltype_mapping, by="Identity")

print("Reference with mapping:")
print(reference)

# Process to tidy format
processed.data <- data %>%
    pivot_longer(-Gene, names_to="name", values_to="value") %>%
    filter(!str_detect(name, "[.]"))

# Create plots for each gene
dir.create("plots", showWarnings = FALSE)

for (k in unique(data$Gene)) {
  cat("\nProcessing gene:", k, "\n")
  
  # Filter for one gene and merge with reference
  single.data <- processed.data %>%
    filter(Gene == k) %>%
    merge(reference, by.x="name", by.y="Cluster", all.x=TRUE)
  
  cat("Expression range:", min(single.data$value, na.rm=TRUE), "to", max(single.data$value, na.rm=TRUE), "\n")
  
  # Merge with ggPlantmap using mapped ROI names
  final.table <- ggPlantmap.merge(
    ggPm.At.leaf.crosssection,
    single.data,
    id.x = "ROI.name",
    id.y = "ROI_mapped"
  )
  
  # Calculate appropriate limits for color scale
  min_val <- min(single.data$value, na.rm=TRUE)
  max_val <- max(single.data$value, na.rm=TRUE)
  
  # Create heatmap with adjusted scale for negative values
  if (min_val < 0) {
    # Use a diverging color scale for data with negative values
    plot <- ggPlantmap.heatmap(final.table, value.quant = value) + 
      scale_fill_gradient2(
        low="blue", 
        mid="white", 
        high="red", 
        midpoint=0,
        limits=c(min_val, max_val)
      ) +
      labs(title=paste0(k)) +
      theme_minimal()
  } else {
    # Use simple gradient for positive-only data
    plot <- ggPlantmap.heatmap(final.table, value.quant = value) + 
      scale_fill_gradient(low="white", high="red", limits=c(0, max_val)) +
      labs(title=paste0(k)) +
      theme_minimal()
  }
  
  print(plot)
  
  # Save plot
  ggsave(plot=plot, paste0("plots/", k, ".png"), dpi=300, width=8, height=8)
  cat("✓ Saved: plots/", k, ".png\n")
}

cat("\n✓ All plots saved in 'plots' directory\n")

# ===============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(ggPlantmap)
library(xml2)
library(dplyr)

p <- path.expand("~/Desktop/UofT Files/Third Year/BCB330/at.h5ad")
p <- normalizePath(p, mustWork = TRUE)

sce <- readH5AD(p)

# ============================================
# BASIC STRUCTURE
# ============================================

# Overall summary
sce

# Dimensions (cells x genes)
dim(sce)

# ============================================
# GENE INFORMATION
# ============================================

# See all gene names
rownames(sce)  # or head(rownames(sce), 20) for first 20

# Check if AT4G29780 is in the dataset
"SRC2" %in% rownames(sce)

# Gene metadata (if any)
rowData(sce)

# ============================================
# CELL INFORMATION & METADATA
# ============================================

# Cell metadata (this is what you need!)
colData(sce)

# See column names in metadata
colnames(colData(sce))

# View first few rows
head(colData(sce))

# Check unique cell types
unique(colData(sce)$cell_type)  # adjust column name if different

# Check what treatment conditions exist
unique(colData(sce)$condition)  # or treatment, or whatever they named it

# ============================================
# EXPRESSION DATA
# ============================================

# Access the expression matrix
assays(sce)  # see what assays are available

# Get expression data (usually in "X" or "counts" or "logcounts")
expr_matrix <- assay(sce, "X")  # or assay(sce, 1) for first assay

# Dimensions
dim(expr_matrix)

# ============================================
# EXTRACT SRC2 EXPRESSION DATA
# ============================================

# Get expression for SRC2 across all cells
src2_expr <- assay(sce, "X")["SRC2", ]

# Combine with cell metadata using the correct column
gene_data <- data.frame(
  label_v2 = colData(sce)$label_v2,      # detailed cell type
  label_major = colData(sce)$label_major, # major cell type (simpler)
  condition = colData(sce)$condition,
  expression = src2_expr
)

head(gene_data)

# ============================================
# CHECK WHAT CELL TYPES ARE AVAILABLE
# ============================================

# See all unique cell types (detailed)
unique(gene_data$label_v2)

# See all unique major cell types (simpler)
table(gene_data$label_major)

# ============================================
# CALCULATE MEAN EXPRESSION
# ============================================

library(dplyr)

# Mean expression by major cell type (across all conditions)
mean_by_major <- gene_data %>%
  group_by(label_major) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE))

print(mean_by_major)

# Mean expression by detailed cell type
mean_by_detailed <- gene_data %>%
  group_by(label_v2) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE))

print(mean_by_detailed)

# Mean by cell type AND condition
mean_by_type_condition <- gene_data %>%
  group_by(label_major, condition) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

print(mean_by_type_condition)

# Create a mapping to broader categories
gene_data_categorized <- gene_data %>%
  mutate(broad_category = case_when(
    grepl("Epidermal", label_v2) ~ "epidermis",
    grepl("Guard", label_v2) ~ "epidermis",  # Guard cells are specialized epidermal
    grepl("Vascular", label_v2) ~ "xylem",   # Assuming vascular includes xylem
    grepl("Phloem|Sieve element", label_v2) ~ "phloem",
    TRUE ~ "other"
  ))

# Calculate mean expression for your three target cell types
mean_for_heatmap <- gene_data_categorized %>%
  filter(broad_category %in% c("epidermis", "xylem", "phloem")) %>%
  group_by(broad_category) %>%
  summarize(mean_expr = mean(expression, na.rm = TRUE))

print(mean_for_heatmap)

# Create external_data for your heatmap
external_data <- tibble(
  ROI.name = mean_for_heatmap$broad_category,
  Gene.expression = exp(mean_for_heatmap$mean_expr)
)

print(external_data)

# Now merge and create heatmap
merged_data <- ggPlantmap.merge(
  vascular_data,
  external_data,
  id.x = "ROI.name",
  id.y = "ROI.name"
)

ggPlantmap.heatmap(merged_data, Gene.expression) +
  scale_fill_gradient(low = "white", high = "red")

# ================================================

# Function to convert Icy XML to ggPlantmap format
# Function to convert Icy XML to ggPlantmap format
icy_to_ggplantmap <- function(xml_file) {
  
  # Read XML file
  xml_data <- read_xml(xml_file)
  
  # Find all ROI nodes
  roi_nodes <- xml_find_all(xml_data, ".//rois/roi")
  
  if (length(roi_nodes) == 0) {
    stop("No ROI nodes found in the XML file. Check the file structure.")
  }
  
  # Parse each ROI
  roi_list <- lapply(seq_along(roi_nodes), function(roi_index) {
    roi <- roi_nodes[[roi_index]]
    
    # Get ROI name from the <name> tag (biological names like C1, C2, END, etc.)
    roi_name <- xml_text(xml_find_first(roi, ".//name"))
    if (is.na(roi_name) || roi_name == "") {
      roi_name <- paste0("ROI_", roi_index)
    }
    
    # Get all points
    point_nodes <- xml_find_all(roi, ".//points/point")
    
    if (length(point_nodes) == 0) {
      return(NULL)
    }
    
    # Extract coordinates for each point
    coords <- data.frame(
      ROI.name = roi_name,
      ROI.id = roi_index,
      point = seq_along(point_nodes),
      x = sapply(point_nodes, function(p) {
        as.numeric(xml_text(xml_find_first(p, ".//pos_x")))
      }),
      y = sapply(point_nodes, function(p) {
        as.numeric(xml_text(xml_find_first(p, ".//pos_y")))
      })
    )
    
    return(coords)
  })
  
  # Remove NULL entries and combine
  roi_list <- roi_list[!sapply(roi_list, is.null)]
  
  if (length(roi_list) == 0) {
    stop("No valid ROI data found in the XML file.")
  }
  
  # Combine all ROIs into single data frame
  roi_data <- bind_rows(roi_list)
  
  # Convert to tibble for better printing
  roi_data <- as_tibble(roi_data)
  
  return(roi_data)
}

# Use it on your file
roi_data <- icy_to_ggplantmap("/Users/stevenqiao/Desktop/ggplantmap_data/icy_new_3.xml")

# View the result
print(roi_data)

ggPlantmap.plot(roi_data)

# =======================================================
# Sample external quantitative data for heatmap testing
# This should match your cell type names from your SVG
library(tibble)

vascular_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/vascular.xml")

ggPlantmap.to.SVG(vascular_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "vascular_output.svg")

guard_cell_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/guard_cell.xml")

ggPlantmap.to.SVG(guard_cell_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "guard_cell_output.svg")

epidermal_cell_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/epidermal_cell.xml")

ggPlantmap.to.SVG(epidermal_cell_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "epidermal_cell_output.svg")

trichome_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/trichome.xml")

ggPlantmap.to.SVG(trichome_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "trichome_output.svg")

palisade_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/palisade.xml")

ggPlantmap.to.SVG(palisade_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "palisade_output.svg")

spongy_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/spongy.xml")

ggPlantmap.to.SVG(spongy_data,
                  group.name = "ROI.name",
                  author = "Steven Qiao",
                  svg.name = "spongy_output.svg")



# Generate random gene expression values for all 34 regions
external_data <- tibble(
  ROI.name = c("epidermis", "xylem", "phloem"),
  Gene.expression = c(100, 0, 50)  # Random values between 0-100
)

print(external_data)

merged_data <- ggPlantmap.merge(
  vascular_data, 
  external_data,
  id.x = "ROI.name",
  id.y = "ROI.name"
)

library(ggplot2)

print(merged_data)

ggPlantmap.heatmap(merged_data, Gene.expression) +
  scale_fill_gradient(low = "white", high = "red")

# =======================================================

vascular_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/vascular.xml")
print(vascular_data)
ggPlantmap.plot(vascular_data)

epidermal_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/epidermal.xml")
print(epidermal_data)
ggPlantmap.plot(epidermal_data)

pavement_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/pavement.xml")
print(pavement_data)
ggPlantmap.plot(pavement_data)

trichome_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/trichome.xml")
print(trichome_data)
ggPlantmap.plot(trichome_data)

palisade_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/palisade.xml")
print(palisade_data)
ggPlantmap.plot(palisade_data)

spongy_data <- XML.to.ggPlantmap("/Users/stevenqiao/Desktop/spongy.xml")
print(spongy_data)
ggPlantmap.plot(spongy_data)


library(ggplot2)
library(dplyr)

# ============================================
# EXTRACT UMAP COORDINATES
# ============================================

# Get UMAP coordinates from reducedDims
umap_coords <- reducedDim(sce, "X_umap")

# Combine UMAP coordinates with cell metadata and SRC2 expression
umap_data <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  cell_type_detailed = colData(sce)$label_v2,
  cell_type_major = colData(sce)$label_major,
  condition = colData(sce)$condition,
  SRC2_expression = assay(sce, "X")["SRC2", ]
)

# ============================================
# PLOT 1: UMAP BY CELL TYPE (DETAILED)
# ============================================

ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cell_type_detailed)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "UMAP colored by Cell Type",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

# ============================================
# PLOT 2: UMAP BY MAJOR CELL TYPE (SIMPLER)
# ============================================

ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cell_type_major)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "UMAP colored by Major Cell Type",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Major Cell Type"
  )

# ============================================
# PLOT 3: UMAP BY SRC2 EXPRESSION
# ============================================

ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = SRC2_expression)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_gradient2(
    low = "blue",
    mid = "lightgray",
    high = "red",
    midpoint = median(umap_data$SRC2_expression, na.rm = TRUE)
  ) +
  theme_minimal() +
  labs(
    title = "UMAP colored by SRC2 Expression",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "SRC2 Expression"
  )

# ============================================
# ALTERNATIVE: VIRIDIS COLOR SCALE FOR SRC2
# ============================================

ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = SRC2_expression)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    title = "UMAP colored by SRC2 Expression (Viridis)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "SRC2 Expression"
  )








