import xml.etree.ElementTree as ET


def convert_to_eplant_format(input_svg, output_svg):
    # Register the SVG namespace without a prefix
    ET.register_namespace('', 'http://www.w3.org/2000/svg')

    # Parse the SVG
    tree = ET.parse(input_svg)
    root = tree.getroot()

    # Define namespace
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    # Find all groups
    for group in root.findall('.//svg:g', ns):
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
            path.set('id', polygon.get('id'))

            # Apply group attributes to path
            for attr, value in group_attrs.items():
                path.set(attr, value)

            # Add path to group and remove polygon
            group.append(path)
            group.remove(polygon)

        # Remove styling attributes from group (keep only id)
        for attr in list(group_attrs.keys()):
            if attr in group.attrib:
                del group.attrib[attr]

    # Write output
    tree.write(output_svg, encoding='utf-8', xml_declaration=True)


# Use it
convert_to_eplant_format('spongy_output.svg', 'spongy_eplant_format.svg')
