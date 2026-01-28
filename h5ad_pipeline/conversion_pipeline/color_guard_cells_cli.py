#!/usr/bin/env python3
"""
Color guard cells in the merged SVG based on expression data.
This script processes D0, W0, R15, and W15 guard cells.

Usage:
    python color_guard_cells_cli.py <svg_file> <json_file> <output_file> [options]

Example:
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg --min -0.2 --max 0.2
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg --opacity 0.9
"""

import json
import argparse
import sys
from xml.etree import ElementTree as ET


def expression_to_color(value, min_val=-0.5, max_val=0.5):
    """
    Convert expression value to color using a yellow to dark red gradient.

    For log-normalized data (which is often negative):
    - Higher values (less negative, closer to 0) = HIGH expression = DARK RED
    - Lower values (more negative, further from 0) = LOW expression = YELLOW

    This provides high contrast between values while keeping everything visible.

    Args:
        value: Expression value (can be negative for log-normalized data)
        min_val: Minimum expression value (maps to yellow - lowest expression)
        max_val: Maximum expression value (maps to dark red - highest expression)

    Returns:
        Hex color string
    """
    # Normalize value to 0-1 range
    if value <= min_val:
        normalized = 0.0  # Lowest expression -> yellow
    elif value >= max_val:
        normalized = 1.0  # Highest expression -> dark red
    else:
        # CRITICAL: This correctly maps higher values to 1.0 (dark red)
        normalized = (value - min_val) / (max_val - min_val)

    # Yellow (#FFD700) to Orange (#FF8C00) to Red (#FF0000) to Dark Red (#8B0000)
    # This creates a heat map effect with high contrast

    if normalized < 0.33:
        # Yellow to Orange
        t = normalized / 0.33
        r = 255
        g = int(215 - (215 - 140) * t)  # 215 (gold) to 140 (orange)
        b = int(0)
    elif normalized < 0.67:
        # Orange to Red
        t = (normalized - 0.33) / 0.34
        r = 255
        g = int(140 * (1 - t))  # 140 (orange) to 0 (red)
        b = 0
    else:
        # Red to Dark Red
        t = (normalized - 0.67) / 0.33
        r = int(255 - (255 - 139) * t)  # 255 (red) to 139 (dark red)
        g = 0
        b = 0

    return f"#{r:02x}{g:02x}{b:02x}"


def load_expression_data(json_file):
    """Load expression data from JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Create a dictionary mapping cell_type to expression value
    expression_dict = {}
    for entry in data:
        cell_type = entry['cell_type']
        expression_dict[cell_type] = entry['average_expression']

    return expression_dict


def color_guard_cells(svg_file, expression_dict, output_file, min_val=-0.5, max_val=0.5, opacity=0.7):
    """
    Color guard cells in SVG based on expression data.

    Args:
        svg_file: Input SVG file path
        expression_dict: Dictionary mapping cell_type to expression values
        output_file: Output SVG file path
        min_val: Minimum value for color scale (lowest expression)
        max_val: Maximum value for color scale (highest expression)
        opacity: Fill opacity (0.0 to 1.0)
    """
    # Parse SVG
    tree = ET.parse(svg_file)
    root = tree.getroot()

    # Define namespace
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    # Target conditions and their guard cell types
    conditions = ['D0', 'W0', 'R15', 'W15']

    # Collect all guard cell expression values first to show proper context
    guard_values = {}
    for condition in conditions:
        cell_type = f"{condition}_Guard"
        if cell_type in expression_dict:
            guard_values[cell_type] = expression_dict[cell_type]

    if not guard_values:
        print("Warning: No guard cell expression data found!")
        return

    # Determine actual min/max from data if not specified
    actual_min = min(guard_values.values())
    actual_max = max(guard_values.values())

    print("=" * 70)
    print("Guard Cell Expression Values and Colors")
    print("=" * 70)
    print(f"Data range: {actual_min:.4f} to {actual_max:.4f}")
    print(f"Color scale: {min_val:.4f} (YELLOW=low) to {max_val:.4f} (DARK RED=high)")
    print(f"Opacity: {opacity}")
    print("-" * 70)
    print(f"{'Cell Type':<20} {'Expression':<12} {'Color':<10} {'Interpretation'}")
    print("-" * 70)

    # Sort by expression value (highest to lowest) to show ranking
    sorted_cells = sorted(guard_values.items(), key=lambda x: x[1], reverse=True)

    for cell_type, expression_value in sorted_cells:
        color = expression_to_color(expression_value, min_val, max_val)

        # Determine interpretation
        if expression_value == actual_max:
            interp = "HIGHEST expression"
        elif expression_value == actual_min:
            interp = "LOWEST expression"
        else:
            interp = "Medium expression"

        print(f"{cell_type:<20} {expression_value:>8.4f}    {color:<10} {interp}")

    print("=" * 70)

    # Now apply colors to SVG
    colored_count = 0
    for condition in conditions:
        cell_type = f"{condition}_Guard"

        # Get expression value
        if cell_type not in expression_dict:
            continue

        expression_value = expression_dict[cell_type]
        color = expression_to_color(expression_value, min_val, max_val)

        # Find guard cell group by ID
        guard_group_id = f"{condition}_guard"

        # Find all path elements within guard cell groups
        for elem in root.iter():
            # Check if this is a group with the guard ID
            if elem.get('id') == guard_group_id:
                # Find all path elements within this group
                for path in elem.iter('{http://www.w3.org/2000/svg}path'):
                    # Set fill color instead of none
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    colored_count += 1

    # Write output
    tree.write(output_file, encoding='utf-8', xml_declaration=True)
    print(f"\n✓ Colored {colored_count} paths across {len(guard_values)} guard cell types")
    print(f"✓ Output written to: {output_file}")
    print()


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Color guard cells in SVG based on expression data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic usage (auto-detect data range):
    %(prog)s merged_grid.svg avg_expression.json output.svg

  With custom color scale:
    %(prog)s merged_grid.svg avg_expression.json output.svg --min -0.2 --max 0.2

  With custom opacity:
    %(prog)s merged_grid.svg avg_expression.json output.svg --opacity 0.9

  Combined options:
    %(prog)s merged_grid.svg avg_expression.json output.svg --min -0.15 --max 0.15 --opacity 0.8

Note:
  For log-normalized data, higher values (less negative) indicate higher expression.
  The color scale uses a heat map: yellow (low) → orange → red → dark red (high).
  This provides high contrast between different expression levels.
        """
    )

    # Required arguments
    parser.add_argument('svg_file', help='Input SVG file path')
    parser.add_argument('json_file', help='Input JSON file with expression data')
    parser.add_argument('output_file', help='Output SVG file path')

    # Optional arguments
    parser.add_argument('--min', type=float, default=None,
                        help='Minimum expression value for color scale (default: auto from data)')
    parser.add_argument('--max', type=float, default=None,
                        help='Maximum expression value for color scale (default: auto from data)')
    parser.add_argument('--opacity', type=float, default=0.7,
                        help='Fill opacity from 0.0 to 1.0 (default: 0.7)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress output messages')

    # Parse arguments
    args = parser.parse_args()

    # Validate opacity
    if not 0.0 <= args.opacity <= 1.0:
        print("Error: Opacity must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)

    try:
        # Load expression data
        if not args.quiet:
            print("Loading expression data...")
        expression_dict = load_expression_data(args.json_file)

        # Auto-detect min/max from guard cell data if not specified
        guard_values = []
        for condition in ['D0', 'W0', 'R15', 'W15']:
            cell_type = f"{condition}_Guard"
            if cell_type in expression_dict:
                guard_values.append(expression_dict[cell_type])

        if not guard_values:
            print("Error: No guard cell expression data found in JSON", file=sys.stderr)
            sys.exit(1)

        # Set min/max from data if not specified
        min_val = args.min if args.min is not None else min(guard_values)
        max_val = args.max if args.max is not None else max(guard_values)

        # Validate min/max
        if min_val >= max_val:
            print("Error: --min must be less than --max", file=sys.stderr)
            sys.exit(1)

        # Color guard cells
        if not args.quiet:
            print("Coloring guard cells...\n")
        color_guard_cells(args.svg_file, expression_dict, args.output_file,
                          min_val, max_val, args.opacity)

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
