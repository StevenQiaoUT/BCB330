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

1. **Convert ICY XML to SVG format:**
```r
# Run in R
source("icy_xml_to_ggplantmap_svg.R")
```

2. **Process H5AD data:**
```bash
python umap.py
```

3. **Generate composite grid layout:**
```bash
python fourby3tree.py
```

## Workflow

### Step 1: Prepare ICY XML Files
Place your ICY XML cell annotation files in the root directory. Each file should represent a specific cell type (epidermal, guard, trichome, palisade, spongy, vascular).

### Step 2: Convert XML to SVG
```r
# Open R and run the conversion script
source("icy_xml_to_ggplantmap_svg.R")
```

This generates individual SVG files for each cell type and timepoint.

### Step 3: Process Single-Cell Data
```bash
# Generate UMAP visualizations
python umap.py

# View H5AD data (optional)
Rscript h5ad_viewer_ad.r
```

### Step 4: Create Composite Visualization
```bash
# Generate 4×3 grid layout (4 timepoints × 3 tissue categories)
python fourby3tree.py
```

This creates `merged_grid.svg` with all cell types organized by tissue category across timepoints.

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
- **vascular**: Vascular bundle (includes BS, xylem, phloem)

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
