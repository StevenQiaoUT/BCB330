#!/usr/bin/env python3
"""
SVG to XML Coordinate Extractor
Reads an XML file with empty coords attributes and fills them based on SVG layout.
"""

import re
import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom


def calculate_bbox(base_x_min, base_y_min, base_x_max, base_y_max, translate_x, translate_y, scale=0.3):
    """
    Calculate bounding box after applying transform.

    Args:
        base_x_min, base_y_min, base_x_max, base_y_max: Original bounding box coordinates
        translate_x, translate_y: Translation offsets
        scale: Scale factor (default 0.3)

    Returns:
        String in format: "minx,miny,maxx,miny,maxx,maxy,minx,maxy"
    """
    # Apply scale first, then translate
    x_min = base_x_min * scale + translate_x
    y_min = base_y_min * scale + translate_y
    x_max = base_x_max * scale + translate_x
    y_max = base_y_max * scale + translate_y

    # Format as rectangle corners: minx,miny,maxx,miny,maxx,maxy,minx,maxy
    return f"{x_min:.0f},{y_min:.0f},{x_max:.0f},{y_min:.0f},{x_max:.0f},{y_max:.0f},{x_min:.0f},{y_max:.0f}"


def generate_coordinates():
    """
    Generate all coordinates based on SVG layout.

    Returns:
        Dictionary mapping tissue_id to coordinate string
    """
    # Define base bounding boxes for each cell type (from analyzing the SVG paths)
    base_boxes = {
        'bsc': (98, 42, 429, 398),  # Bundle sheath cells
        'phloem': (172, 206, 349, 315),  # Phloem
        'xylem': (167, 107, 352, 198),  # Xylem
        'guard_cell': (21, 17, 86, 97),  # Guard cells
        'epidermal_cell': (10, 4, 86, 96),  # Epidermal cells
        'trichome': (13, 12, 88, 100),  # Trichomes
        'spongy': (8, 6, 100, 79),  # Spongy mesophyll
        'palisade': (13, 9, 80, 82)  # Palisade mesophyll
    }

    # Position offsets (translate values for columns)
    positions = {
        'D0': 0,  # Moderate drought
        'R15': 200,  # Recovered
        'W0': 400,  # Well watered
        'W15': 600  # Well watered and irrigation
    }

    # Row offsets for each cell type
    row_offsets = {
        'bsc': 0,  # Vascular tissues at top
        'phloem': 0,  # Vascular tissues at top
        'xylem': 0,  # Vascular tissues at top
        'guard_cell': 200,  # Row 2
        'epidermal_cell': 400,  # Row 3
        'trichome': 600,  # Row 4
        'spongy': 800,  # Row 5
        'palisade': 1000  # Row 6
    }

    # Calculate coordinates for all combinations
    coordinates = {}

    for cell_type, (base_x_min, base_y_min, base_x_max, base_y_max) in base_boxes.items():
        for position_name, col_offset in positions.items():
            # Calculate translation
            translate_x = col_offset
            translate_y = row_offsets[cell_type]

            # Generate coordinates
            coords = calculate_bbox(base_x_min, base_y_min, base_x_max, base_y_max,
                                    translate_x, translate_y)

            tissue_id = f"{cell_type}_{position_name}"
            coordinates[tissue_id] = coords

    return coordinates


def prettify_xml(elem):
    """
    Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, encoding='utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t", encoding='utf-8').decode('utf-8')


def fill_coordinates(input_file, output_file):
    """
    Read XML file, fill in empty coords attributes, and write to output file.

    Args:
        input_file: Path to input XML file with empty coords
        output_file: Path to output XML file with filled coords
    """
    # Generate all coordinates
    coordinates = generate_coordinates()

    # Parse the input XML file
    tree = ET.parse(input_file)
    root = tree.getroot()

    # Find all tissue elements and update their coords
    tissues_updated = 0
    tissues_not_found = []

    for tissue in root.findall('.//tissue'):
        tissue_id = tissue.get('id')

        if tissue_id in coordinates:
            # Find the area element and update coords
            area = tissue.find('area')
            if area is not None:
                area.set('coords', coordinates[tissue_id])
                tissues_updated += 1
                print(f"Updated {tissue_id}: {coordinates[tissue_id]}")
            else:
                print(f"Warning: No <area> element found for {tissue_id}")
        else:
            tissues_not_found.append(tissue_id)
            print(f"Warning: No coordinates generated for {tissue_id}")

    # Write the updated XML to output file
    # First convert to string with pretty printing
    xml_string = prettify_xml(root)

    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(xml_string)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"Summary:")
    print(f"  Tissues updated: {tissues_updated}")
    if tissues_not_found:
        print(f"  Tissues not found: {len(tissues_not_found)}")
        for tid in tissues_not_found:
            print(f"    - {tid}")
    print(f"  Output written to: {output_file}")
    print(f"{'=' * 60}")


def main():

    input_file = "merged_grid.xml"
    output_file = "merged_grid_coords.xml"

    try:
        fill_coordinates(input_file, output_file)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except ET.ParseError as e:
        print(f"Error: Failed to parse XML file: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
