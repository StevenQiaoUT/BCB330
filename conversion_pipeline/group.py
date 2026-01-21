import xml.etree.ElementTree as ET


def restructure_svg_by_tissue(input_svg, output_svg):
    # Register namespace
    ET.register_namespace('', 'http://www.w3.org/2000/svg')

    # Parse the SVG
    tree = ET.parse(input_svg)
    root = tree.getroot()

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

    # Find all existing groups and categorize their paths
    groups_to_remove = []

    for group in root.findall('.//svg:g', ns):
        group_id = group.get('id')

        if group_id:
            # Determine tissue type
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

                # Mark group for removal
                groups_to_remove.append(group)

    # Remove old groups
    for group in groups_to_remove:
        root.remove(group)

    # Add new tissue groups to root
    for tissue_group in tissue_groups.values():
        if len(tissue_group) > 0:  # Only add if it has paths
            root.append(tissue_group)

    # Write output
    tree.write(output_svg, encoding='utf-8', xml_declaration=True)


# Use it
restructure_svg_by_tissue('spongy_eplant_format.svg', 'spongy_tissue_grouped.svg')
