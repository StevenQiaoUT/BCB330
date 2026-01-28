"""
Composite Tissue Grid Creator

Create a 4×3 grid where each cell contains multiple cell types grouped by tissue category.

Layout:
- Row 1 (Epidermal): Guard Cell, Epidermal cell, Trichome
- Row 2 (Mesophyll): Palisade, Spongy
- Row 3 (Vascular): Composite circular structure

Usage:
    python fourbythree.py <input_dir> <output_file> [--cols COLS] [--cell-width WIDTH] [--cell-height HEIGHT] [--spacing SPACING]

Arguments:
    input_dir    : Directory containing input SVG files with naming pattern {timepoint}_{celltype}.svg
    output_file  : Path to output merged SVG file

Options:
    --cols         : Number of columns (timepoints) (default: 4)
    --cell-width   : Width of each grid cell in pixels (default: 400)
    --cell-height  : Height of each grid cell in pixels (default: 300)
    --spacing      : Space between cells in pixels (default: 50)

Examples:
    # Basic usage with default settings
    python fourbythree.py conditioned/ merged_grid.svg

    # Custom dimensions
    python fourbythree.py conditioned/ output.svg --cell-width 500 --cell-height 400

    # Full custom configuration
    python fourbythree.py data/ grid.svg --cols 4 --cell-width 400 --cell-height 300 --spacing 50
"""

import sys
import os
import re
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path


def create_composite_tissue_grid(svg_files_dict, output_file, cols=4, cell_width=400, cell_height=300, spacing=50):
    """
    Create a 4×3 grid where each cell contains multiple cell types grouped by tissue category.

    Parameters:
    - svg_files_dict: nested dict {timepoint: {cell_type: filepath}}
    - output_file: path to output SVG file
    - cols: number of columns (timepoints)
    - cell_width/cell_height: dimensions of each grid cell
    - spacing: space between cells
    """

    # Define tissue organization
    tissue_layout = {
        'epidermal': ['guard', 'epidermal', 'trichome'],
        'mesophyll': ['palisade', 'spongy'],
        'vascular': ['vascular']
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
            'vascular': {'x': 100, 'scale': 0.4}
        }
    }

    timepoints = sorted(svg_files_dict.keys())
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

            # Add a background rectangle for visualization
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
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(svg_content))
        print(f'✓ Created 4×3 tissue grid in {output_file}')
        print(f'  Layout: {rows} tissue categories × {cols} timepoints')
        print(f'  Total dimensions: {total_width} × {total_height}')
    except Exception as e:
        print(f'Error writing output file: {e}')
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Create a composite tissue grid from SVG files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('input_dir', help='Directory containing input SVG files')
    parser.add_argument('output_file', help='Path to output merged SVG file')
    parser.add_argument('--cols', type=int, default=4, help='Number of columns (default: 4)')
    parser.add_argument('--cell-width', type=int, default=400, help='Width of each grid cell (default: 400)')
    parser.add_argument('--cell-height', type=int, default=300, help='Height of each grid cell (default: 300)')
    parser.add_argument('--spacing', type=int, default=50, help='Space between cells (default: 50)')

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f'Error: Input directory "{args.input_dir}" does not exist')
        sys.exit(1)

    # Ensure input_dir ends with /
    input_dir = args.input_dir if args.input_dir.endswith('/') else args.input_dir + '/'

    print("=" * 60)
    print("Composite Tissue Grid Creator")
    print("=" * 60)
    print(f"Input directory:  {input_dir}")
    print(f"Output file:      {args.output_file}")
    print(f"Grid layout:      3 rows × {args.cols} columns")
    print(f"Cell dimensions:  {args.cell_width} × {args.cell_height}")
    print(f"Cell spacing:     {args.spacing}")
    print("=" * 60)
    print()

    # Define expected timepoints and cell types
    timepoints = ['D0', 'R15', 'W0', 'W15']
    cell_types = ['guard', 'epidermal', 'trichome', 'palisade', 'spongy', 'vascular']

    # Build the file dictionary
    svg_files_dict = {}
    files_found = 0
    files_missing = 0

    for timepoint in timepoints:
        svg_files_dict[timepoint] = {}
        for cell_type in cell_types:
            filename = f'{input_dir}{timepoint}_{cell_type}.svg'
            if Path(filename).exists():
                svg_files_dict[timepoint][cell_type] = filename
                files_found += 1
            else:
                print(f'Warning: {filename} not found')
                files_missing += 1

    print(f'\nFile Discovery:')
    print(f'  Found: {files_found} files')
    if files_missing > 0:
        print(f'  Missing: {files_missing} files')
    print()

    if files_found == 0:
        print('Error: No SVG files found in input directory')
        print(f'Expected file pattern: {{timepoint}}_{{celltype}}.svg')
        print(f'Example: D0_guard.svg, W0_palisade.svg')
        sys.exit(1)

    # Create the composite grid
    create_composite_tissue_grid(
        svg_files_dict,
        args.output_file,
        cols=args.cols,
        cell_width=args.cell_width,
        cell_height=args.cell_height,
        spacing=args.spacing
    )

    print('\nDone!')


if __name__ == "__main__":
    main()
