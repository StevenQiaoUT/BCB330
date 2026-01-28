# BCB330 - Single-Cell RNA-seq Visualization Pipeline

Pipeline for converting ICY XML format to ggPlantmap-compatible SVG format and processing H5AD single-cell RNA-seq data for visualization in ePlant.

## Overview

This project processes single-cell RNA-seq data from Arabidopsis leaf tissues across different drought stress timepoints (W0, D0, R15, W15) and generates interactive SVG visualizations for integration into the ePlant browser. The pipeline handles multiple cell types including guard cells, epidermal cells, trichomes, palisade mesophyll, spongy mesophyll, and vascular tissues.

## Table of Contents
- [Features](#features)
- [Installation Instructions](#installation-instructions)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Workflow](#workflow)
- [File Naming Conventions](#file-naming-conventions)
- [Output](#output)

## Features

- Convert ICY XML cell annotations to ggPlantmap-compatible SVG format
- Process H5AD single-cell RNA-seq data files
- Generate UMAP visualizations for different cell types
- Create composite tissue grid layouts for multiple timepoints
- Support for multiple Arabidopsis leaf tissue types

## Installation Instructions

This project requires Python, R, and system-level dependencies. Follow these steps in order.


### Prerequisites

- **Python**: 3.8 or higher
- **R**: 4.5.x (tested with R 4.5.1)
- **Git**: For cloning repositories


### Step 1: Install System Dependencies

#### HDF5 Library

HDF5 is required for reading `.h5ad` files and must be installed before Python packages.

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install libhdf5-dev
```

**macOS (with Homebrew):**
```bash
brew install hdf5
```

**Windows:**
```bash
# Using Conda (recommended for Windows)
conda install -c conda-forge hdf5
```

### Step 2: Install Python Dependencies

```bash
pip install -r requirements.txt
```

### Step 3: Install R Dependencies

#### 3.1 Update BiocManager and Bioconductor

**IMPORTANT**: You must update BiocManager to the latest version and install Bioconductor 3.22 for R 4.5 compatibility.

Open R or RStudio and run:

```r
install.packages("BiocManager")

BiocManager::install(version = "3.22")
```

#### 3.2 Install Bioconductor Packages

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("zellkonverter")
```

#### 3.3 Install tidyverse

```r
install.packages("tidyverse")
```

#### 3.4 Install ggPlantmap from GitHub

```r
install.packages("devtools")
library(devtools)

install_github("leonardojo/ggPlantmap")
```

### Version Information

This setup has been tested with:
- R 4.5.1
- Bioconductor 3.22
- Python 3.8+

## Project Structure
```
BCB330/
├── **h5ad_pipeline/**
│   ├── h5ad_to_json.py
│   ├── h5ad_viewer_ad.r
│   ├── h5ad_viewer_bc.r
│   ├── umap.py
│   ├── **conversion_pipeline/**
│   │   ├── fourbythree.py
│   │   ├── get_coordinates.py
│   │   ├── gg_to_eplant_ungrouped.py
│   │   ├── group.py
│   │   ├── combine.py
│   |── **svg_generation/**
│   │   ├── icy_xml_to_ggplantmap_svg.r
│   │   ├── merged_grid.svg
│   │   ├── merged_grid.xml
│   │   ├── merged_grid_coords.xml
│   │   ├── eplant_config_with_coords.xml
│   │   ├── conditioned/
│   │   ├── icy_outputs/
│   │   └── intermediates/
├── .Rproj.user/
├── .gitignore
├── BCB330 Proposal - Steven Qiao.pdf
├── ePlant SVG and Expression Data Guide.pdf
├── README.md
├── requirements.txt
└── timer.r
```

## Usage

### Quick Start

#### Generate UMAP Plots
```bash
# Basic usage - visualize a specific gene
python interactive_umap.py at.h5ad AT3G05727

# Custom output filename
python interactive_umap.py at.h5ad AT3G05727 my_umap.html

# Use files from different directory
python interactive_umap.py /path/to/data.h5ad AT3G05727 output.html
```

#### Convert H5AD Expression Values to JSON
```bash
# Single gene
python h5ad_to_json_avg.py at.h5ad output.json label_majorXcondition AT3G05727

# Multiple genes
python h5ad_to_json_avg.py at.h5ad output.json label_majorXcondition AT3G05727 AT1G01010 AT5G12345

# All genes (no gene list specified)
python h5ad_to_json_avg.py at.h5ad all_genes.json label_v2

# Different paths
python h5ad_to_json_avg.py /data/at.h5ad /output/result.json cell_type AT3G05727
```

### Pipeline Usage

#### 1. Trace PNG/JPG Using ICY
Visit: https://icy.bioimageanalysis.org/

#### 2. Convert ICY XML to ggPlantmap SVG format
```r
# Open R and run the conversion script
source("icy_xml_to_ggplantmap_svg.r")

# Convert one XML file to SVG
Rscript icy_xml_to_ggplantmap_svg.r vascular.xml vascular_output.svg

# With author name
Rscript icy_xml_to_ggplantmap_svg.r guard_cell.xml guard_output.svg --author "Steven Qiao"

# Using full paths
Rscript icy_xml_to_ggplantmap_svg.r /path/to/input/epidermal.xml /path/to/output/epidermal.svg
```

#### 3. Convert from ggPlantmap SVG to ePlant SVG Format
```bash
# Basic usage
python gg_to_eplant_ungrouped.py spongy_output.svg spongy_eplant_format.svg

# Process file from different directory
python gg_to_eplant_ungrouped.py data/input.svg output/result.svg
```

#### 4. Group each ePlant SVG by Cell Subtype
```bash
# Basic usage
python group.py spongy_eplant_format.svg spongy_tissue_grouped.svg

# Different paths
python group.py data/input.svg output/grouped.svg
```

#### 5. Combine the Individual 24 ePlant SVGs into a 4×3 Grid
```bash
# Basic usage with default settings
python fourbythree.py conditioned/ merged_grid.svg

# Custom dimensions
python fourbythree.py conditioned/ output.svg --cell-width 500 --cell-height 400

# Full custom configuration
python fourbythree.py data/ grid.svg --cols 4 --cell-width 400 --cell-height 300 --spacing 50
```

#### 6. Generate Coordinates for Grid (Optional)
```bash
# Basic usage
python get_coordinates.py merged_grid.xml merged_grid_coords.xml

# Custom scale factor
python get_coordinates.py input.xml output.xml --scale 0.5
```

## Workflow

### Step 1: Prepare ICY XML Files
Place your ICY XML cell annotation files in the root directory after tracing each one from its original PNG/JPG format. Each file should represent a specific cell type (epidermal, guard, trichome, palisade, spongy, vascular).

In this repository, the ICY XMLs reside in ```conversion_pipeline/icy_outputs```

### Step 2: Convert ICY XML to ggPlantmap SVG
```r
# Open R and run the conversion script
source("icy_xml_to_ggplantmap_svg.r")
```

This generates individual SVG files for each cell type. The process involves first converting the ICY XML into ggPlantmap. Afterwards, ggPlantmap gets converted to ggPlantmap SVG
- Program Input: ```conversion_pipeline/icy_outputs/*.xml```
- Program Output: ```conversion_pipeline/intermediates/*_output.svg```

### Step 3: Convert from ggPlantmap SVG to ePlant SVG Format
```bash
python gg_to_eplant_ungrouped.py spongy_output.svg spongy_eplant_format.svg
```
This step is responsible for rearranging the SVG drawing structure so that it is compatible with ePlant.
- Program Input: ```conversion_pipeline/intermediates/*_output.svg```
- Program Output: ```conversion_pipeline/intermediates/*_eplant_format.svg```

### Step 4: Group by Subtype within Major Cell Type
```bash
python group.py spongy_eplant_format.svg spongy_tissue_grouped.svg
```
For the vascular cell type, the SVG is grouped into bundle sheath cells, xylem, and phloem. The rest of the cells are homogenous in terms of cell type
- Program Input: ```conversion_pipeline/intermediates/*_eplant_format.svg```
- Program Output: ```conversion_pipeline/conditioned/*.svg``` (4 copies for each condition and 6 total cell types yield 24 individual SVGs)

### Step 5: 4×3 Grid Formation
```bash
python fourbythree.py conditioned/ merged_grid.svg
```
The 24 individual SVGs are combined to form a 4×3 grid with condition as column and cell type as rows
- Program Input: ```conversion_pipeline/conditioned/*.svg```
- Program Output: ```conversion_pipeline/merged_grid.svg```

## File Naming Conventions

### Timepoints
- **W0**: Well-watered control
- **D0**: Drought stress initiation
- **R15**: Recovery
- **W15**: Well-watered

### Cell Types
- **guard**: Guard cells
- **epidermal**: Epidermal cells
- **trichome**: Trichome cells
- **palisade**: Palisade mesophyll
- **spongy**: Spongy mesophyll
- **vascular**: Vascular bundle (includes bundle sheath cells, xylem, phloem)

### File Format
Individual SVG files follow the pattern: `{timepoint}_{cell_type}.svg`

Example: `W0_guard.svg`, `D0_palisade.svg`

## Output

### Individual SVG Files
Located in `conditioned/` directory, one file per cell type per timepoint.

### Composite Grid
- **File**: `merged_grid.svg`
- **Layout**: 4 columns (W0, D0, R15, W15) × 3 rows (Epidermal, Mesophyll, Vascular)
- **Dimensions**: Configurable in `fourbythree.py`

### Visualization Features
- Consistent stroke styling (1px, dark gray)
- Scaled cell representations for size uniformity
- Labeled tissue categories and timepoints
- Grid structure for easy comparison across conditions

### Interactive UMAP Visualizations
- **File**: `interactive_umap.html`
- HTML-based interactive plots with dropdown menus for cell type selection
- Color-coded gene expression levels
- Hover tooltips with cell type and expression information
- Statistics for each cell type (mean, standard deviation, cell count)

### JSON Expression Data
- **File**: `avg_expression.json`
- Average gene expression values per cell type
- Compatible with ePlant browser integration
- Supports single or multiple genes

## Troubleshooting

### Common Issues

**Issue**: "cannot change working directory" error
- **Solution**: All scripts now use command-line arguments. Run them from any directory with full paths to input/output files.

**Issue**: SVG files not found
- **Solution**: Ensure ICY XML files are converted first and specify correct input directory path.

**Issue**: Python module not found
- **Solution**: Install missing packages: `pip install [package-name]`

**Issue**: R package errors
- **Solution**: Update BiocManager and install Bioconductor 3.22 (see Installation Instructions)

**Issue**: "packages are not available for this version of R"
- **Solution**: Make sure BiocManager is updated: `install.packages("BiocManager")` then `BiocManager::install(version = "3.22")`

**Issue**: HDF5 library errors
- **Solution**: Install HDF5 system library before Python packages (see Step 1 of Installation)

**Issue**: Gene not found in dataset
- **Solution**: Check available genes in your H5AD file. Scripts will suggest available genes if specified gene is not found.

## References

See ggPlantMap ICY guide (https://github.com/leonardojo/ggPlantmap/blob/main/guides/TutorialforXMLfile.pdf) for detailed specifications on ICY XML formation.

## License

MIT License

## Contact

Steven Qiao - University of Toronto
- Project: BCB330Y (Bioinformatics and Computational Biology)
- Lab: Provart Lab
