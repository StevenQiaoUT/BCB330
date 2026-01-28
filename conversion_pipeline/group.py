"""
SVG Tissue Grouping Tool

Restructures SVG by grouping paths by tissue type.
Consolidates multiple tissue-specific groups into single parent groups.

Usage:
    python group.py <input_svg> <output_svg>

Arguments:
    input_svg   : Path to input SVG file with individual cell groups
    output_svg  : Path to output SVG file with tissue-level grouping

Examples:
    # Basic usage
    python group.py spongy_eplant_format.svg spongy_tissue_grouped.svg

    # Process file from different directory
    python group.py data/input.svg output/grouped.svg
"""

import sys
import os
import xml.etree.ElementTree as ET


def restructure_svg_by_tissue(input_svg, output_svg):
    """
    Restructure SVG by grouping paths by tissue type.

    Groups paths with IDs starting with tissue names into tissue-level groups:
    - BSC, phloem, xylem, guard_cell, epidermal_cell, trichome, palisade, spongy

    Args:
        input_svg: Path to input SVG file
        output_svg: Path to output SVG file
    """
    print(f"Loading SVG from: {input_svg}")

    # Register namespace
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

    # Create new tissue groups
    tissue_groups = {
        'BSC': ET.Element('g', {'id': 'BSC'}),
        'phloem': ET.Element('g', {'id': 'phloem'}),
        'xylem': ET.Element('g', {'id': 'xylem'}),
        'guard_cell': ET.Element('g', {'id': 'guard_cell'}),
        'epidermal_cell': ET.Element('g', {'id': 'epidermal_cell'}),
        'trichome': ET.Element('g', {'id': 'trichome'}),
        'palisade': ET.Element('g', {'id': 'palisade'}),
        'spongy': ET.Element('g', {'id': 'spongy'})
    }

    # Track statistics
    tissue_counts = {tissue: 0 for tissue in tissue_groups.keys()}
    groups_to_remove = []

    # Find all existing groups and categorize their paths
    for group in root.findall('.//svg:g', ns):
        group_id = group.get('id')

        if group_id:
            # Determine tissue type based on group ID prefix
            tissue_type = None
            if group_id.startswith('BSC'):
                tissue_type = 'BSC'
            elif group_id.startswith('phloem'):
                tissue_type = 'phloem'
            elif group_id.startswith('xylem'):
                tissue_type = 'xylem'
            elif group_id.startswith('guard_cell'):
                tissue_type = 'guard_cell'
            elif group_id.startswith('epidermal_cell'):
                tissue_type = 'epidermal_cell'
            elif group_id.startswith('trichome'):
                tissue_type = 'trichome'
            elif group_id.startswith('palisade'):
                tissue_type = 'palisade'
            elif group_id.startswith('spongy'):
                tissue_type = 'spongy'

            if tissue_type:
                # Move all paths from this group to the tissue group
                for path in group.findall('svg:path', ns):
                    # Keep the individual cell ID on the path
                    if not path.get('id'):
                        path.set('id', group_id)
                    tissue_groups[tissue_type].append(path)
                    tissue_counts[tissue_type] += 1

                # Mark group for removal
                groups_to_remove.append(group)

    # Remove old groups
    for group in groups_to_remove:
        root.remove(group)

    # Add new tissue groups to root
    tissues_added = 0
    for tissue_name, tissue_group in tissue_groups.items():
        if len(tissue_group) > 0:  # Only add if it has paths
            root.append(tissue_group)
            tissues_added += 1

    # Write output
    print(f"Writing output to: {output_svg}")
    try:
        tree.write(output_svg, encoding='utf-8', xml_declaration=True)
        print(f"\n✓ Restructuring successful!")
        print(f"  Groups removed:      {len(groups_to_remove)}")
        print(f"  Tissue groups added: {tissues_added}")
        print(f"\nPaths per tissue:")
        for tissue_name, count in tissue_counts.items():
            if count > 0:
                print(f"  {tissue_name:20s}: {count:4d} paths")
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
    print("SVG Tissue Grouping Tool")
    print("=" * 60)
    print(f"Input:  {input_svg}")
    print(f"Output: {output_svg}")
    print("=" * 60)
    print()

    restructure_svg_by_tissue(input_svg, output_svg)

    print("\nDone!")


if __name__ == "__main__":
    main()
