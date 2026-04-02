"""
GGPlantmap to ePlant Format Converter

Converts SVG polygons to paths and moves group attributes to child elements
for ePlant compatibility.

Usage:
    python gg_to_eplant_ungrouped.py <input_svg> <output_svg>

Arguments:
    input_svg   : Path to input SVG file (from ggPlantmap)
    output_svg  : Path to output SVG file (ePlant format)

Examples:
    # Basic usage
    python gg_to_eplant_ungrouped.py spongy_output.svg spongy_eplant_format.svg

    # Process file from different directory
    python gg_to_eplant_ungrouped.py data/input.svg output/result.svg
"""

import sys
import os
import xml.etree.ElementTree as ET


def convert_to_eplant_format(input_svg, output_svg):
    """
    Convert ggPlantmap SVG to ePlant format.

    - Converts polygons to paths
    - Moves group styling attributes to child elements
    - Maintains element IDs

    Args:
        input_svg: Path to input SVG file
        output_svg: Path to output SVG file
    """
    print(f"Loading SVG from: {input_svg}")

    # Register the SVG namespace without a prefix
    ET.register_namespace('', 'http://www.w3.org/2000/svg')

    # Parse the SVG
    try:
        tree = ET.parse(input_svg)
        root = tree.getroot()
    except Exception as e:
        print(f"Error parsing SVG file: {e}")
        sys.exit(1)

    # Define namespace
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    polygons_converted = 0
    groups_processed = 0

    # Find all groups
    for group in root.findall('.//svg:g', ns):
        groups_processed += 1

        # Get group attributes that should be moved to children
        group_attrs = {
            'fill': group.get('fill'),
            'stroke': group.get('stroke'),
            'stroke-width': group.get('stroke-width'),
            'stroke-linecap': group.get('stroke-linecap'),
            'stroke-linejoin': group.get('stroke-linejoin'),
            'stroke-opacity': group.get('stroke-opacity')
        }

        # Remove None values
        group_attrs = {k: v for k, v in group_attrs.items() if v is not None}

        # Process each polygon in the group
        for polygon in list(group.findall('svg:polygon', ns)):
            # Get points
            points = polygon.get('points')

            if not points:
                print(f"  Warning: Polygon without points attribute in group {group.get('id')}")
                continue

            # Convert polygon points to path d attribute
            # Points format: "x1,y1 x2,y2 x3,y3..."
            points_list = points.strip().split()
            path_d = 'M' + points_list[0].replace(',', ' ')
            for point in points_list[1:]:
                path_d += ' L' + point.replace(',', ' ')
            path_d += ' Z'  # Close the path

            # Create new path element (without namespace prefix)
            path = ET.Element('path')
            path.set('d', path_d)

            # Preserve polygon ID if it exists
            polygon_id = polygon.get('id')
            if polygon_id:
                path.set('id', polygon_id)

            # Apply group attributes to path
            for attr, value in group_attrs.items():
                path.set(attr, value)

            # Add path to group and remove polygon
            group.append(path)
            group.remove(polygon)
            polygons_converted += 1

        # Remove styling attributes from group (keep only id and transform)
        for attr in list(group_attrs.keys()):
            if attr in group.attrib:
                del group.attrib[attr]

    # Write output
    print(f"Writing output to: {output_svg}")
    try:
        tree.write(output_svg, encoding='utf-8', xml_declaration=True)
        print(f"\n✓ Conversion successful!")
        print(f"  Groups processed:    {groups_processed}")
        print(f"  Polygons converted:  {polygons_converted}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_svg = sys.argv[1]
    output_svg = sys.argv[2]

    # Validate input file
    if not os.path.exists(input_svg):
        print(f"Error: Input file '{input_svg}' not found")
        sys.exit(1)

    if not input_svg.lower().endswith('.svg'):
        print(f"Warning: Input file '{input_svg}' does not have .svg extension")

    print("=" * 60)
    print("GGPlantmap to ePlant Format Converter")
    print("=" * 60)
    print(f"Input:  {input_svg}")
    print(f"Output: {output_svg}")
    print("=" * 60)
    print()

    convert_to_eplant_format(input_svg, output_svg)

    print("\nDone!")


if __name__ == "__main__":
    main()
