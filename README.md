# BCB330 - Single-Cell RNA-seq Visualization Pipeline

Pipeline for converting ICY XML format to ggPlantmap-compatible SVG format and processing H5AD single-cell RNA-seq data for visualization in ePlant.

## Overview

This project processes single-cell RNA-seq data from Arabidopsis leaf tissues across different drought stress timepoints (W0, D0, R15, W15) and generates interactive SVG visualizations for integration into the ePlant browser. The pipeline handles multiple cell types including guard cells, epidermal cells, trichomes, palisade mesophyll, spongy mesophyll, and vascular tissues.

## Table of Contents
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
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

## Prerequisites

### Software Requirements
- **R** (4.0+)
- **Python** (3.8+)

### R Packages
```r
install.packages(c("ggplot2", "xml2", "dplyr", "ggPlantmap", "SingleCellExperiment", "zellkonverter"))
```

### Python Packages
```bash
pip install -r requirements.txt
```

Required Python packages likely include:
- pandas
- numpy
- scanpy
- plotly
- pathlib

## Installation
```bash
# Clone the repository
git clone https://github.com/StevenQiaoUT/BCB330.git
cd BCB330

# Install Python dependencies
pip install -r requirements.txt

# Install R packages (run in R console)
# See Prerequisites section above
```

## Project Structure
```
BCB330/
├── conversion_pipeline/
│   ├── conditioned/
│   ├── icy_outputs/                 
│   ├── intermediates/
│   ├── combine.py                 
│   ├── icy_xml_to_ggplantmap_svg.r  
│   ├── fourbythree.py        
|   ├── get_coordinates.py        
│   ├── gg_to_eplant_ungrouped.py
│   └── group.py    
├── h5ad_viewer_*.r
├── h5ad_to_json.py          
├── umap.py                   
├── timer.r
├── requirements.txt                 
├── BCB330 Proposal - Steven Qiao.pdf
└── ePlant SVG and Expression Data Guide.pdf
```

## Usage

### Quick Start
#### Generate UMAP Plots
```
umap.py
```
#### Convert h5ad Expressions Values to JSON
```
h5ad_to_json.py
```

#### Pipeline Usage
1. **Trace PNG/JPG Using ICY:**
```
https://icy.bioimageanalysis.org/
```

2. **Convert ICY XML to ggPlantmap SVG format:**
```
icy_xml_to_ggplantmap_svg.r
```

2. **Convert from ggPlantmap SVG to ePlant SVG Format:**
```
gg_to_ePlant_ungrouped.py
```

3. **Group each ePlant SVG by Cell Subtype:**
```
group.py
```

4. **Combine the Individual 24 ePlant SVGs into a 4*3 Grid**
```
fourbythree.py
```

## Workflow

### Step 1: Prepare ICY XML Files
Place your ICY XML cell annotation files in the root directory after tracing each one from its original PNG/JPG format. Each file should represent a specific cell type (epidermal, guard, trichome, palisade, spongy, vascular).

In this repository, the ICY XMLs reside in ```conversion_pipeline/icy_outputs```

### Step 2: Convert ICY XML to ggPlantmap SVG
```
# Open R and run the conversion script
source("icy_xml_to_ggplantmap_svg.r")
```

This generates individual SVG files for each cell type. The process involves first converting the ICY XML into ggPlantmap. Afterwards, ggPlantmap gets converted to ggPlantmap SVG
- Program Input: conversion_pipeline/icy_outputs/*.xml
- Program Output: conversion_pipeline/intermediates/*_output.svg

### Step 3: Convert from ggPlantmap SVG to ePlant SVG Format
```
gg_to_ePlant_ungrouped.py
```
This step is responsible for rearranging the SVG drawing structure so that it is compatible with ePlant.
- Program Input: conversion_pipeline/intermediates/*_output.svg
- Program Output: conversion_pipeline/intermediates/*_eplant_format.svg

### Step 4: Group by Subtype within Major Cell Type
```
group.py
```
For the vascular cell type, the SVG is grouped into bundle sheath cells, xylem, and phloem. The rest of the cells are homogenous in terms of cell type
- Program Input: conversion_pipeline/intermediates/*_eplant_format.svg
- Program Output: conversion_pipeline/conditioned/*.svg (4 copies for each condition and 6 total cell types yield 24 individual SVGs)

### Step 5: 4x3 Grid Formation
```
group.py
```
The 24 individual SVGs are combined to form a 4x3 grid with condition as column and cell type as rows
- Program Input: conversion_pipeline/conditioned/*.svg
- Program Output: conversion_pipeline/merged_grid.svg

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

## Troubleshooting

### Common Issues

**Issue**: SVG files not found
- **Solution**: Ensure ICY XML files are converted first and output files are in `conditioned/` directory

**Issue**: Python module not found
- **Solution**: Install missing packages: `pip install [package-name]`

**Issue**: R package errors
- **Solution**: Install required R packages (see Prerequisites)

## References

See ggPlantMap ICY guide (https://github.com/leonardojo/ggPlantmap/blob/main/guides/TutorialforXMLfile.pdf) for detailed specifications on ICY XML formation.

## License

MIT License

## Contact

Steven Qiao - University of Toronto
- Project: BCB330Y (Bioinformatics and Computational Biology)
- Lab: Provart Lab
