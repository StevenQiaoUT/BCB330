import xml.etree.ElementTree as ET
from pathlib import Path
import re


def create_svg_grid(svg_files, output_file, cols=6, spacing=200, scale_factor=None):
    """
    Merge SVG files into a grid layout while preserving group IDs.

    Parameters:
    - scale_factor: dict mapping tissue types to scale values, or single float for all
    """

    # Calculate grid dimensions
    rows = (len(svg_files) + cols - 1) // cols

    # Get dimensions from first SVG
    first_svg = ET.parse(svg_files[0])
    first_root = first_svg.getroot()

    # Try to get viewBox or default dimensions
    viewbox = first_root.get('viewBox', '0 0 600 600').split()
    svg_width = float(viewbox[2])
    svg_height = float(viewbox[3])

    # Set total dimensions
    total_width = cols * spacing
    total_height = rows * spacing

    # Start building SVG as string for more control
    svg_content = [
        '<?xml version="1.0" encoding="utf-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_width} {total_height}" width="{total_width}" height="{total_height}">'
    ]

    # Process each SVG file
    for idx, svg_file in enumerate(svg_files):
        # Calculate grid position
        row = idx // cols
        col = idx % cols
        x_offset = col * spacing
        y_offset = row * spacing

        # Read and parse SVG content
        with open(svg_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Extract everything between <svg> and </svg>
        # Remove XML declaration and svg tags
        content = re.sub(r'<\?xml.*?\?>', '', content)
        content = re.sub(r'<svg[^>]*>', '', content)
        content = re.sub(r'</svg>', '', content)

        # Remove any xmlns attributes from child elements
        content = re.sub(r'\s+xmlns="[^"]*"', '', content)

        # Determine scale for this tissue type
        tissue_type = Path(svg_file).stem.split('_')[1]  # Extract tissue name

        if scale_factor is not None:
            if isinstance(scale_factor, dict):
                scale = scale_factor.get(tissue_type, 1.0)
            else:
                scale = scale_factor
        else:
            scale = 1.0

        # Build transform with both translate and scale
        if scale != 1.0:
            transform = f'translate({x_offset}, {y_offset}) scale({scale})'
        else:
            transform = f'translate({x_offset}, {y_offset})'

        # Wrap in a positioned group
        svg_content.append(f'  <g id="position_{Path(svg_file).stem}" transform="{transform}">')
        svg_content.append(content.strip())
        svg_content.append('  </g>')

    svg_content.append('</svg>')

    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(svg_content))

    print(f'✓ Merged {len(svg_files)} SVGs into {output_file}')
    print(f'  Grid layout: {rows} rows × {cols} columns')
    print(f'  Total dimensions: {total_width} × {total_height}')


# Define your SVG files
tissue_types = ['vascular', 'guard', 'epidermal', 'trichome', 'spongy', 'palisade']
timepoints = ['D0', 'R15', 'W0', 'W15']

# Option 1: Tissues as rows, timepoints as columns (like your example)
svg_files = []
for tissue in tissue_types:
    for timepoint in timepoints:
        filename = f'conditioned/{timepoint}_{tissue}.svg'
        if Path(filename).exists():
            svg_files.append(filename)
        else:
            print(f'Warning: {filename} not found')

print(f'Found {len(svg_files)} SVG files')

# Define scale factors for each tissue type
# Adjust these values to make all cells similar size
tissue_scales = {
    'vascular': 0.3,  # Scale down vascular cells significantly
    'guard': 1.0,  # Keep guard cells at normal size
    'epidermal': 1.0,  # Keep epidermal at normal size
    'trichome': 1.0,  # Keep trichome at normal size
    'spongy': 1.0,  # Keep spongy at normal size
    'palisade': 1.0  # Keep palisade at normal size
}

# Create grid: 4 columns (timepoints), 6 rows (tissues)
create_svg_grid(svg_files, 'merged_grid.svg', cols=4, spacing=200, scale_factor=tissue_scales)
