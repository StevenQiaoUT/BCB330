import xml.etree.ElementTree as ET
from pathlib import Path
import re


def create_composite_tissue_grid(svg_files_dict, output_file, cols=4, cell_width=400, cell_height=300, spacing=50):
    """
    Create a 4×3 grid where each cell contains multiple cell types grouped by tissue category.

    Layout:
    - Row 1 (Epidermal): Guard Cell, Epidermal cell, Trichome (left to right)
    - Row 2 (Mesophyll): Palisade, Spongy (left to right)
    - Row 3 (Vascular): Composite circular structure (outer ring, BS, xylem, phloem)

    - Columns: 4 timepoints (W0, D0, R15, W15)

    Parameters:
    - svg_files_dict: nested dict {timepoint: {cell_type: filepath}}
    - cols: number of columns (timepoints)
    - cell_width/cell_height: dimensions of each grid cell
    - spacing: space between cells
    """

    # Define tissue organization
    tissue_layout = {
        'epidermal': ['guard', 'epidermal', 'trichome'],
        'mesophyll': ['palisade', 'spongy'],
        'vascular': ['vascular']  # This one contains multiple subtypes in one file
    }

    rows = 3  # epidermal, mesophyll, vascular

    # Calculate total dimensions
    total_width = cols * (cell_width + spacing) - spacing
    total_height = rows * (cell_height + spacing) - spacing

    # Start building SVG
    svg_content = [
        '<?xml version="1.0" encoding="utf-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_width} {total_height}" '
        f'width="{total_width}" height="{total_height}">'
    ]

    # Define positions for cell types within each tissue category
    # These control horizontal spacing of cell types within one grid cell
    cell_type_positions = {
        'epidermal': {
            'guard': {'x': 0, 'scale': 1.2},
            'epidermal': {'x': 130, 'scale': 1.2},
            'trichome': {'x': 260, 'scale': 1.2}
        },
        'mesophyll': {
            'palisade': {'x': 50, 'scale': 1.3},
            'spongy': {'x': 220, 'scale': 1.3}
        },
        'vascular': {
            'vascular': {'x': 100, 'scale': 0.4}  # Scaled down from 0.6
        }
    }

    timepoints = ['W0', 'D0', 'R15', 'W15']
    tissue_categories = ['epidermal', 'mesophyll', 'vascular']

    # Process each grid cell
    for row_idx, tissue_category in enumerate(tissue_categories):
        for col_idx, timepoint in enumerate(timepoints):

            # Calculate grid cell position
            grid_x = col_idx * (cell_width + spacing)
            grid_y = row_idx * (cell_height + spacing)

            # Create a group for this grid cell
            svg_content.append(
                f'  <g id="{tissue_category}_{timepoint}" '
                f'transform="translate({grid_x}, {grid_y})">'
            )

            # Add a background rectangle for visualization (optional, can remove)
            svg_content.append(
                f'    <rect x="0" y="0" width="{cell_width}" height="{cell_height}" '
                f'fill="none" stroke="#4a4a4a" stroke-width="1" stroke-dasharray="5,5"/>'
            )

            # Add label for this cell
            svg_content.append(
                f'    <text x="{cell_width / 2}" y="20" text-anchor="middle" '
                f'font-family="Arial" font-size="14" font-weight="bold">'
                f'{tissue_category.capitalize()} - {timepoint}</text>'
            )

            # Add each cell type for this tissue category
            cell_types = tissue_layout[tissue_category]

            for cell_type in cell_types:
                # Get the SVG file for this timepoint and cell type
                if timepoint in svg_files_dict and cell_type in svg_files_dict[timepoint]:
                    svg_file = svg_files_dict[timepoint][cell_type]

                    if Path(svg_file).exists():
                        # Get positioning info
                        pos_info = cell_type_positions[tissue_category][cell_type]
                        cell_x = pos_info['x']
                        cell_scale = pos_info['scale']
                        cell_y = 50  # Vertical offset from top of cell

                        # Read SVG content
                        with open(svg_file, 'r', encoding='utf-8') as f:
                            content = f.read()

                        # Clean up content
                        content = re.sub(r'<\?xml.*?\?>', '', content)
                        content = re.sub(r'<svg[^>]*>', '', content)
                        content = re.sub(r'</svg>', '', content)
                        content = re.sub(r'\s+xmlns="[^"]*"', '', content)

                        # Update stroke attributes to dark gray and width 1
                        content = re.sub(r'stroke="[^"]*"', 'stroke="#4a4a4a"', content)
                        content = re.sub(r'stroke-width="[^"]*"', 'stroke-width="1"', content)

                        # Add cell type group with positioning
                        transform = f'translate({cell_x}, {cell_y}) scale({cell_scale})'
                        svg_content.append(
                            f'    <g id="{timepoint}_{cell_type}" transform="{transform}">'
                        )
                        svg_content.append(content.strip())
                        svg_content.append('    </g>')
                    else:
                        print(f'Warning: {svg_file} not found')
                else:
                    print(f'Warning: No file for {timepoint} - {cell_type}')

            # Close grid cell group
            svg_content.append('  </g>')

    svg_content.append('</svg>')

    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(svg_content))

    print(f'✓ Created 4×3 tissue grid in {output_file}')
    print(f'  Layout: {rows} tissue categories × {cols} timepoints')
    print(f'  Total dimensions: {total_width} × {total_height}')


# Build the file dictionary
timepoints = ['D0', 'R15', 'W0', 'W15']
cell_types = ['guard', 'epidermal', 'trichome', 'palisade', 'spongy', 'vascular']

svg_files_dict = {}
for timepoint in timepoints:
    svg_files_dict[timepoint] = {}
    for cell_type in cell_types:
        filename = f'conditioned/{timepoint}_{cell_type}.svg'
        if Path(filename).exists():
            svg_files_dict[timepoint][cell_type] = filename
        else:
            print(f'Warning: {filename} not found')

# Create the composite grid
create_composite_tissue_grid(
    svg_files_dict,
    'merged_grid.svg',
    cols=4,  # 4 timepoints
    cell_width=400,  # Width of each grid cell
    cell_height=300,  # Height of each grid cell
    spacing=50  # Space between cells
)
