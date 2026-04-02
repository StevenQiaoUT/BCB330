#!/usr/bin/env Rscript
#' ICY XML to ggPlantmap SVG Converter
#'
#' Converts ICY XML cell annotation files to ggPlantmap-compatible SVG format.
#'
#' Usage:
#'   Rscript icy_xml_to_ggplantmap_svg.r <input_xml> <output_svg> [--author AUTHOR]
#'
#' Arguments:
#'   input_xml   : Path to input ICY XML file
#'   output_svg  : Path to output SVG file
#'   --author    : (Optional) Author name for SVG metadata (default: "User")
#'
#' Examples:
#'   # Basic usage
#'   Rscript icy_xml_to_ggplantmap_svg.r vascular.xml vascular_output.svg
#'
#'   # With author name
#'   Rscript icy_xml_to_ggplantmap_svg.r guard_cell.xml guard_output.svg --author "Steven Qiao"
#'
#'   # Process multiple files with a batch script
#'   for file in *.xml; do
#'     Rscript icy_xml_to_ggplantmap_svg.r "$file" "${file%.xml}_output.svg"
#'   done

# ============================================
# LOAD REQUIRED LIBRARIES
# ============================================

# Suppress package startup messages
suppressPackageStartupMessages({
  library(ggPlantmap)
  library(xml2)
  library(dplyr)
  library(tibble)
})

# ============================================
# PARSE COMMAND-LINE ARGUMENTS
# ============================================

args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
print_usage <- function() {
  cat("\n")
  cat("ICY XML to ggPlantmap SVG Converter\n")
  cat("===================================\n\n")
  cat("Usage:\n")
  cat("  Rscript icy_xml_to_ggplantmap_svg.r <input_xml> <output_svg> [--author AUTHOR]\n\n")
  cat("Arguments:\n")
  cat("  input_xml   : Path to input ICY XML file\n")
  cat("  output_svg  : Path to output SVG file\n")
  cat("  --author    : (Optional) Author name for SVG metadata (default: 'User')\n\n")
  cat("Examples:\n")
  cat("  Rscript icy_xml_to_ggplantmap_svg.r vascular.xml vascular_output.svg\n")
  cat("  Rscript icy_xml_to_ggplantmap_svg.r guard.xml guard.svg --author 'Steven Qiao'\n\n")
}

# Check for minimum arguments
if (length(args) < 2) {
  print_usage()
  stop("Error: Insufficient arguments provided", call. = FALSE)
}

# Parse arguments
input_xml <- args[1]
output_svg <- args[2]
author <- "User"  # default

# Check for optional --author flag
if (length(args) >= 4 && args[3] == "--author") {
  author <- args[4]
}

# Validate input file exists
if (!file.exists(input_xml)) {
  stop(paste("Error: Input file", input_xml, "not found"), call. = FALSE)
}

# Print configuration
cat("\n")
cat("============================================================\n")
cat("ICY XML to ggPlantmap SVG Converter\n")
cat("============================================================\n")
cat("Input file:  ", input_xml, "\n")
cat("Output file: ", output_svg, "\n")
cat("Author:      ", author, "\n")
cat("============================================================\n")
cat("\n")

# ============================================
# CONVERSION FUNCTION
# ============================================

icy_to_ggplantmap <- function(xml_file) {
  #' Convert ICY XML to ggPlantmap format
  #'
  #' @param xml_file Path to ICY XML file
  #' @return tibble with ROI data (ROI.name, ROI.id, point, x, y)

  cat("Reading XML file...\n")

  # Read XML file
  tryCatch({
    xml_data <- read_xml(xml_file)
  }, error = function(e) {
    stop(paste("Error reading XML file:", e$message), call. = FALSE)
  })

  # Find all ROI nodes
  roi_nodes <- xml_find_all(xml_data, ".//rois/roi")

  if (length(roi_nodes) == 0) {
    stop("No ROI nodes found in the XML file. Check the file structure.", call. = FALSE)
  }

  cat("Found", length(roi_nodes), "ROI(s)\n")

  # Parse each ROI
  roi_list <- lapply(seq_along(roi_nodes), function(roi_index) {
    roi <- roi_nodes[[roi_index]]

    # Get ROI name from the <name> tag
    roi_name <- xml_text(xml_find_first(roi, ".//name"))
    if (is.na(roi_name) || roi_name == "") {
      roi_name <- paste0("ROI_", roi_index)
    }

    # Get all points
    point_nodes <- xml_find_all(roi, ".//points/point")

    if (length(point_nodes) == 0) {
      warning(paste("ROI", roi_index, "has no points, skipping"))
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
    stop("No valid ROI data found in the XML file.", call. = FALSE)
  }

  # Combine all ROIs into single data frame
  roi_data <- bind_rows(roi_list)

  # Convert to tibble
  roi_data <- as_tibble(roi_data)

  cat("Extracted", nrow(roi_data), "points from", length(unique(roi_data$ROI.name)), "unique ROI(s)\n")

  return(roi_data)
}

# ============================================
# MAIN CONVERSION
# ============================================

cat("Converting ICY XML to ggPlantmap format...\n")

# Convert XML to ggPlantmap format
tryCatch({
  roi_data <- icy_to_ggplantmap(input_xml)
}, error = function(e) {
  stop(paste("Conversion failed:", e$message), call. = FALSE)
})

# Display summary
cat("\nROI Summary:\n")
roi_summary <- roi_data %>%
  group_by(ROI.name) %>%
  summarize(
    points = n(),
    x_range = paste(round(min(x), 1), "-", round(max(x), 1)),
    y_range = paste(round(min(y), 1), "-", round(max(y), 1)),
    .groups = "drop"
  )
print(roi_summary)

# ============================================
# CONVERT TO SVG
# ============================================

cat("\nGenerating SVG...\n")

tryCatch({
  ggPlantmap.to.SVG(
    roi_data,
    group.name = "ROI.name",
    author = author,
    svg.name = output_svg
  )

  cat("\n✓ Successfully saved SVG to:", output_svg, "\n")

  # Verify output file was created
  if (file.exists(output_svg)) {
    file_size <- file.info(output_svg)$size
    cat("  File size:", format(file_size, big.mark = ","), "bytes\n")
  } else {
    warning("Output file was not created successfully")
  }

}, error = function(e) {
  stop(paste("SVG generation failed:", e$message), call. = FALSE)
})

cat("\n============================================================\n")
cat("Conversion complete!\n")
cat("============================================================\n\n")
