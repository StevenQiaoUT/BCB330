# BCB330 - A Pseudo-Bulk Tissue eFP Browser for Single Cell Data for ePlant

A pipeline for visualizing single-cell RNA-seq data from *Arabidopsis thaliana* leaf tissues within the ePlant browser. This project bridges single-cell and bulk tissue expression data by generating interactive SVG-based eFP visualizations and UMAP plots for integration into the BAR ePlant platform.

## Overview

This project processes single-cell RNA-seq data across four drought stress conditions (W0, D0, R15, W15) and six cell types (guard, epidermal, trichome, palisade, spongy, vascular). The pipeline converts ICY XML cell annotations into ePlant-compatible SVG diagrams, assembles them into a 4×3 condition-by-cell-type grid, and serves the result alongside an interactive UMAP via a Python CGI application hosted on the BAR server.

## Table of Contents
- [Features](#features)
- [Installation Instructions](#installation-instructions)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Workflow](#workflow)
- [File Naming Conventions](#file-naming-conventions)
- [Output](#output)
- [Troubleshooting](#troubleshooting)

## Features

- Convert ICY XML cell annotations to ggPlantmap-compatible SVG format
- Transform ggPlantmap SVGs into ePlant-compatible format
- Color SVG cell regions by average gene expression per cell type
- Assemble 24 individual SVGs into a 4×3 grid (condition × cell type)
- Generate interactive Plotly UMAP visualizations with gene expression overlay
- Serve the combined eFP + UMAP viewer via Apache CGI on the BAR server

## Installation Instructions

### Prerequisites

- **Python**: 3.12 or higher (developed with 3.12.1)
- **R**: 4.5.x (tested with R 4.5.1)
- **Git**: For cloning repositories

### Step 1: Install System Dependencies

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
conda install -c conda-forge hdf5
```

### Step 2: Install Python Dependencies

```bash
pip install -r requirements.txt
```

### Step 3: Install R Dependencies

#### 3.1 Update BiocManager and Bioconductor

Open R or RStudio and run:

```r
install.packages("BiocManager")
BiocManager::install(version = "3.22")
```

#### 3.2 Install Bioconductor Packages

```r
BiocManager::install("SingleCellExperiment")
BiocManager::install("zellkonverter")
```

#### 3.3 Install tidyverse and ggPlantmap

```r
install.packages("tidyverse")
install.packages("devtools")
library(devtools)
install_github("leonardojo/ggPlantmap")
```

### Version Information

- R 4.5.1
- Bioconductor 3.22
- Python 3.12+

## Project Structure

```
BCB330/
├── BCB330 Proposal - Steven Qiao.pdf
├── README.md
├── app
│   └── mini_app.py                        # Main CGI viewer (eFP SVG + UMAP)
├── ePlant SVG and Expression Data Guide.pdf
├── pipeline
│   ├── 1_svg_generation                   # ICY XML → ggPlantmap → ePlant SVGs
│   │   ├── icy_xml_to_ggplantmap_svg.r
│   │   ├── icy_outputs/                   # Raw ICY XML files (one per cell type)
│   │   ├── intermediates/                 # Per-cell-type intermediate SVGs
│   │   └── conditioned/                   # Final 24 SVGs (4 conditions × 6 cell types)
│   ├── 2_coordinate_extraction            # Extract SVG coordinates for ePlant XML config
│   │   ├── get_coordinates.py
│   │   └── merged_grid_coords.xml
│   ├── 3_grid_assembly                    # Combine 24 SVGs into 4×3 grid
│   │   ├── color.py
│   │   ├── combine.py
│   │   ├── fourbythree.py
│   │   ├── gg_to_eplant_ungrouped.py
│   │   ├── group.py
│   │   ├── merged_grid.xml
│   │   ├── new_grid.svg
│   │   └── output.svg
│   ├── 4_h5ad_processing                  # scRNA-seq data exploration and export
│   │   ├── h5ad_to_json.py
│   │   ├── h5ad_viewer_ad.r
│   │   ├── h5ad_viewer_bc.r
│   │   └── umap.py
│   └── eplant_config_with_coords.xml      # Final ePlant XML configuration
├── requirements.txt
└── timer.r
```

## Usage

### Running the eFP + UMAP Viewer (BAR Server)

The main application is `main_app/mini_app.py`, deployed as a CGI script on the BAR development server:

```
http://142.150.215.219/~sqiao/cgi-bin/mini_app.py?file=/mnt/home/sqiao/Natanella_integrated_object.h5ad&col=label_majorXcondition&gene=AT3G05727&svg=svg_template.svg
```

**Query Parameters:**

| Parameter | Description | Example |
|-----------|-------------|---------|
| `file` | Absolute path to `.h5ad` file | `/mnt/home/sqiao/Natanella_integrated_object.h5ad` |
| `col` | Cell type metadata column | `label_majorXcondition` |
| `gene` | Arabidopsis gene ID | `AT3G05727` |
| `svg` | SVG template filename | `svg_template.svg` |

**Deploy to bardev:**
```bash
scp -P 12345 -i ~/.ssh/id_ed25519 app/mini_app.py sqiao@142.150.215.219:/mnt/home/sqiao/public_html/cgi-bin/
ssh -p 12345 sqiao@142.150.215.219 "chmod 755 /mnt/home/sqiao/public_html/cgi-bin/mini_app.py"
```

### H5AD Processing

```bash
# Export average expression per cell type to JSON
python pipeline/4_h5ad_processing/h5ad_to_json.py at.h5ad output.json label_majorXcondition AT3G05727

# Multiple genes
python pipeline/4_h5ad_processing/h5ad_to_json.py at.h5ad output.json label_majorXcondition AT3G05727 AT1G01010

# All genes
python pipeline/4_h5ad_processing/h5ad_to_json.py at.h5ad all_genes.json label_v2
```

### Grid Assembly

```bash
# Assemble 24 SVGs into 4×3 grid
python pipeline/3_grid_assembly/fourbythree.py pipeline/1_svg_generation/conditioned/ merged_grid.svg

# Custom dimensions
python pipeline/3_grid_assembly/fourbythree.py conditioned/ output.svg --cell-width 500 --cell-height 400
```

## Workflow

### Step 1: Trace Cell Types Using ICY

Visit [https://icy.bioimageanalysis.org/](https://icy.bioimageanalysis.org/) and trace each cell type from a microscopy image. Export as XML. Place outputs in `pipeline/1_svg_generation/icy_outputs/`.

### Step 2: ICY XML → ggPlantmap SVG

```r
Rscript pipeline/1_svg_generation/icy_xml_to_ggplantmap_svg.r vascular.xml vascular_output.svg
```

- **Input**: `pipeline/1_svg_generation/icy_outputs/*.xml`
- **Output**: `pipeline/1_svg_generation/intermediates/*_output.svg`

### Step 3: ggPlantmap SVG → ePlant SVG Format

```bash
python pipeline/3_grid_assembly/gg_to_eplant_ungrouped.py \
  pipeline/1_svg_generation/intermediates/spongy_output.svg \
  pipeline/1_svg_generation/intermediates/spongy_eplant_format.svg
```

- **Input**: `pipeline/1_svg_generation/intermediates/*_output.svg`
- **Output**: `pipeline/1_svg_generation/intermediates/*_eplant_format.svg`

### Step 4: Group by Cell Subtype

```bash
python pipeline/3_grid_assembly/group.py \
  pipeline/1_svg_generation/intermediates/spongy_eplant_format.svg \
  pipeline/1_svg_generation/conditioned/W0_spongy.svg
```

For vascular tissue, this step groups cells into bundle sheath, xylem, and phloem subtypes. For other cell types the cells are homogeneous. Four copies are made per cell type — one per condition — yielding 24 total SVGs.

- **Input**: `pipeline/1_svg_generation/intermediates/*_eplant_format.svg`
- **Output**: `pipeline/1_svg_generation/conditioned/{condition}_{celltype}.svg`

### Step 5: Assemble 4×3 Grid

```bash
python pipeline/3_grid_assembly/fourbythree.py \
  pipeline/1_svg_generation/conditioned/ \
  pipeline/3_grid_assembly/merged_grid.svg
```

Combines the 24 SVGs into a grid with conditions as columns (W0, D0, R15, W15) and cell types as rows.

- **Input**: `pipeline/1_svg_generation/conditioned/*.svg`
- **Output**: `pipeline/3_grid_assembly/merged_grid.svg`

### Step 6: Extract Coordinates (Optional)

```bash
python pipeline/2_coordinate_extraction/get_coordinates.py \
  pipeline/3_grid_assembly/merged_grid.xml \
  pipeline/2_coordinate_extraction/merged_grid_coords.xml
```

Extracts SVG path coordinates for generating the ePlant XML configuration file.

## File Naming Conventions

### Conditions
| Code | Meaning |
|------|---------|
| `W0` | Well-watered control |
| `D0` | Drought stress |
| `R15` | Drought recovery (15 min after rewatering) |
| `W15` | Well-watered, irrigated 15 min |

### Cell Types
| Code | Cell Type |
|------|-----------|
| `guard` | Guard cells (stomata) |
| `epidermal` | Epidermal cells |
| `trichome` | Trichome cells |
| `palisade` | Palisade mesophyll |
| `spongy` | Spongy mesophyll |
| `vascular` | Vascular bundle (xylem, phloem, bundle sheath) |

### File Format
Individual SVGs follow the pattern: `{condition}_{cell_type}.svg`

Example: `W0_guard.svg`, `R15_palisade.svg`

## Output

### Composite Grid SVG
- **File**: `pipeline/3_grid_assembly/merged_grid.svg`
- **Layout**: 4 columns (W0, D0, R15, W15) × 6 rows (guard, epidermal, trichome, palisade, spongy, vascular)

### Interactive eFP + UMAP Viewer
- Served via Apache CGI at `http://142.150.215.219/~sqiao/cgi-bin/mini_app.py`
- Left panel: ePlant-style SVG colored by average gene expression per cell type
- Right panel: Interactive Plotly UMAP with cell-type dropdown and expression overlay
- Includes gene search bar with Arabidopsis gene ID validation

### JSON Expression Data
- Average gene expression values per cell type, compatible with ePlant browser integration

### ePlant XML Configuration
- **File**: `pipeline/eplant_config_with_coords.xml`
- Maps SVG regions to cell-type labels for ePlant integration

## Troubleshooting

**"cannot change working directory" error**
All scripts use command-line arguments. Run from any directory with full paths to input/output files.

**SVG files not found**
Ensure ICY XML files are converted first and the `conditioned/` directory is populated before running `fourbythree.py`.

**Python module not found**
```bash
pip install -r requirements.txt
```

**R package errors**
Update BiocManager first: `install.packages("BiocManager")` then `BiocManager::install(version = "3.22")`.

**HDF5 library errors**
Install the HDF5 system library before Python packages (see Step 1 of Installation).

**Gene not found in dataset**
Check available genes in your H5AD file. The viewer defaults to `AT3G05727` if no gene is specified.

**Apache CGI timeout**
For large H5AD files, precompute expression averages and cache them. The current deployment uses `Natanella_integrated_object.h5ad` (~144k cells).

## Data Source

Adapted from Dr. Natanella Illouz-Eliaz's publication: "Drought recovery in plants triggers a cell-state-specific immune activation." Dataset: ~144,494 *Arabidopsis thaliana* leaf cells across 16 annotated cell types and 4 drought/recovery conditions, stored in H5AD format with UMAP embeddings.

## References

- ggPlantmap ICY guide: [https://github.com/leonardojo/ggPlantmap/blob/main/guides/TutorialforXMLfile.pdf](https://github.com/leonardojo/ggPlantmap/blob/main/guides/TutorialforXMLfile.pdf)
- ePlant / BAR: [https://bar.utoronto.ca](https://bar.utoronto.ca)
- ICY bioimage analysis: [https://icy.bioimageanalysis.org/](https://icy.bioimageanalysis.org/)

## License

MIT License

## Contact

Steven Qiao — University of Toronto
- Course: BCB330Y (Bioinformatics and Computational Biology Research)
- Lab: Provart Lab, Department of Cell & Systems Biology
