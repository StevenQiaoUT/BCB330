#!/usr/bin/env python3
"""
Color guard cells in the merged SVG based on expression data.
This script processes D0, W0, R15, and W15 guard cells using ePlant's standard color scheme.

Usage:
    python color_guard_cells_cli.py <svg_file> <json_file> <output_file> [options]

Example:
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg --max 2.5
    python color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg --opacity 0.9
"""

import json
import argparse
import sys
from xml.etree import ElementTree as ET


def expression_to_color(value, max_val, scheme='white_to_red'):
    """
    Convert expression value to color.

    Args:
        value: Expression value (already shifted to be non-negative)
        max_val: Maximum expression value for the scale
        scheme: 'white_to_red' (white=low, red=high)
                'yellow_to_red' (yellow=low, red=high)

    Returns:
        Hex color string
    """
    if scheme == 'yellow_to_red':
        minColor = {'red': 255, 'green': 255, 'blue': 0}   # Yellow
    else:  # white_to_red
        minColor = {'red': 255, 'green': 255, 'blue': 255}  # White

    maxColor = {'red': 255, 'green': 0, 'blue': 0}  # Red

    # Calculate ratio relative to maximum
    ratio = value / max_val if max_val > 0 else 0

    # Handle invalid ratios
    if not (0 <= ratio <= float('inf')):
        return "#ffffff"

    # Cap ratio at 1.0
    if ratio > 1:
        ratio = 1

    # Linear interpolation between minColor and maxColor
    red   = minColor['red']   + round((maxColor['red']   - minColor['red'])   * ratio)
    green = minColor['green'] + round((maxColor['green'] - minColor['green']) * ratio)
    blue  = minColor['blue']  + round((maxColor['blue']  - minColor['blue'])  * ratio)

    return f"#{red:02x}{green:02x}{blue:02x}"


def load_expression_data(json_file):
    """Load expression data from JSON file.

    Returns:
        expression_dict: dict mapping cell_type -> average_expression
        gene_name: gene identifier string
        global_max: max single-cell expression from the metadata record (or None)
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    expression_dict = {}
    gene_name = None
    global_max = None
    for entry in data:
        if entry.get('type') == 'metadata':
            gene_name = entry.get('gene', gene_name)
            global_max = entry.get('max_expression', None)
            continue
        cell_type = entry['cell_type']
        expression_dict[cell_type] = entry['average_expression']
        if gene_name is None and 'gene' in entry:
            gene_name = entry['gene']

    return expression_dict, gene_name, global_max


def compute_global_shift(expression_dict):
    """
    Compute the shift needed to make all values non-negative across the
    entire expression_dict.  Returns the shift amount (>= 0).
    """
    all_values = list(expression_dict.values())
    min_val = min(all_values) if all_values else 0
    return -min_val if min_val < 0 else 0


def print_color_table(title, cell_values, global_max, shift, scheme):
    """Print a formatted color table for a set of cell types."""
    actual_min = min(cell_values.values())
    actual_max = max(cell_values.values())
    print("=" * 70)
    print(f"{title}")
    print("=" * 70)
    print(f"Data range (raw):     {actual_min - shift:.4f} to {actual_max - shift:.4f}")
    print(f"Data range (shifted): {actual_min:.4f} to {actual_max:.4f}")
    print(f"Colour ceiling: {global_max:.4f}  (shift applied: {shift:.4f})")
    print(f"Color scheme: {scheme}")
    print("-" * 70)
    print(f"{'Cell Type':<25} {'Shifted val':<14} {'Ratio':<8} {'Color'}")
    print("-" * 70)
    for cell_type, val in sorted(cell_values.items(), key=lambda x: x[1], reverse=True):
        ratio = min(val / global_max, 1.0) if global_max > 0 else 0
        color = expression_to_color(val, global_max, scheme)
        tag = " ← HIGHEST" if val == actual_max else (" ← LOWEST" if val == actual_min else "")
        print(f"{cell_type:<25} {val:>10.4f}    {ratio:>5.3f}    {color}{tag}")
    print("=" * 70)
    print()



def add_legend(root, gene_name, colour_ceiling, shift, n_boxes=10):
    """
    Add an ePlant-style colour scale legend to the bottom-left of the SVG.

    Layout (matches ePlant):
      - Header: gene name (italic) top-left
      - "Single Cell Max" label
      - "Linear" label
      - n_boxes stacked rectangles, red (top) → yellow (bottom)
      - Numeric tick labels on the right of each box boundary
      - Values run from (colour_ceiling - shift) at top to 0 at bottom
        i.e. back-converted to the original unshifted expression space.

    Args:
        root          : root <svg> ET element (modified in-place)
        gene_name     : string used in the header
        colour_ceiling: the shifted max used for coloring (pure red)
        shift         : the shift that was applied to raw values
        n_boxes       : number of colour steps (default 10)
    """
    SVG_NS = "http://www.w3.org/2000/svg"

    # ── Geometry ──────────────────────────────────────────────────────────────
    box_w           = 24
    box_h           = 20
    tick_gap        = 6          # gap between box right edge and tick label
    font_size       = 10
    header_fs       = 11
    left_pad        = 10         # padding left of boxes
    right_pad       = 20         # padding right of tick labels before content starts
    # Approximate width of tick label text (5 chars × ~6px/char)
    tick_label_w    = 40
    legend_width    = left_pad + box_w + tick_gap + tick_label_w + right_pad

    # ── Expand canvas width and shift all existing content rightward ──────────
    SVG_NS_local = "http://www.w3.org/2000/svg"

    # Read current width
    width_attr = root.get("width", "")
    try:
        old_w = float("".join(c for c in width_attr if c.isdigit() or c == "."))
        w_unit = "".join(c for c in width_attr if c.isalpha())
        new_w = old_w + legend_width
        root.set("width", f"{new_w}{w_unit}")
    except ValueError:
        old_w = None

    # Update viewBox width too
    vb = root.get("viewBox", "")
    if vb:
        try:
            vb_parts = vb.split()
            vb_parts[2] = str(float(vb_parts[2]) + legend_width)
            root.set("viewBox", " ".join(vb_parts))
        except (IndexError, ValueError):
            pass

    # Wrap all current children (except any already-added title/legend) in a
    # <g> translated right by legend_width so the legend gets the left margin.
    existing = [c for c in list(root)
                if c.get("id") not in ("expression_legend",)]
    for child in existing:
        root.remove(child)
    shift_g = ET.Element(f"{{{SVG_NS_local}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    # ── Legend position: left margin, bottom-aligned ─────────────────────────
    x0 = left_pad

    # Read canvas height for vertical positioning
    height_attr = root.get("height", "")
    try:
        canvas_h = float("".join(c for c in height_attr if c.isdigit() or c == "."))
    except ValueError:
        try:
            canvas_h = float(root.get("viewBox", "0 0 800 600").split()[3])
        except (IndexError, ValueError):
            canvas_h = 600

    # Total legend height: 3 header lines + n_boxes rows + bottom label
    header_lines = 3          # gene, "Single Cell Max", "Linear"
    header_h     = header_lines * (header_fs + 4)
    legend_h     = header_h + n_boxes * box_h + font_size + 6
    y0           = canvas_h - legend_h - 10   # 10px padding from bottom

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

    # ── Helper: small text element ────────────────────────────────────────────
    def txt(x, y, text, fs=font_size, anchor="start", weight="normal", style=""):
        el = ET.Element(f"{{{SVG_NS}}}text")
        el.set("x", str(x))
        el.set("y", str(y))
        el.set("font-family", "Arial, sans-serif")
        el.set("font-size", str(fs))
        el.set("font-weight", weight)
        el.set("text-anchor", anchor)
        el.set("dominant-baseline", "auto")
        el.set("fill", "#222222")
        if style:
            el.set("font-style", style)
        el.text = text
        return el

    # ── Header labels ─────────────────────────────────────────────────────────
    hy = y0 + header_fs
    g.append(txt(x0, hy, gene_name, fs=header_fs, weight="bold", style="italic"))
    hy += header_fs + 4
    g.append(txt(x0, hy, "Single Cell Max", fs=font_size))
    hy += font_size + 4
    g.append(txt(x0, hy, "Linear", fs=font_size))

    # ── Colour boxes (red at top, yellow at bottom) ───────────────────────────
    boxes_y0 = y0 + header_h

    # The unshifted expression range represented by the colour scale:
    # top    = colour_ceiling - shift  (the raw value that maps to pure red)
    # bottom = 0 - shift               = -shift (the raw value that maps to yellow)
    # But we display 0 at the bottom tick since that is the lowest meaningful value
    # we are showing (yellow = lowest, after shift it is 0).
    raw_top    = colour_ceiling - shift   # unshifted value at pure-red end
    raw_bottom = -shift                   # unshifted value at yellow end (often negative)

    for i in range(n_boxes):
        # i=0 → top box (red), i=n_boxes-1 → bottom box (yellow)
        ratio_top    = 1.0 - i       / n_boxes
        ratio_bottom = 1.0 - (i + 1) / n_boxes
        # Use midpoint ratio for the box fill colour
        ratio_mid = (ratio_top + ratio_bottom) / 2
        shifted_val = ratio_mid * colour_ceiling
        color = expression_to_color(shifted_val, colour_ceiling, "yellow_to_red")

        bx = x0
        by = boxes_y0 + i * box_h

        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x",      str(bx))
        rect.set("y",      str(by))
        rect.set("width",  str(box_w))
        rect.set("height", str(box_h))
        rect.set("fill",   color)
        rect.set("stroke", "#888888")
        rect.set("stroke-width", "0.5")
        g.append(rect)

        # Tick label at the TOP edge of each box (= ratio_top value)
        raw_tick = raw_bottom + ratio_top * (raw_top - raw_bottom)
        tick_label = f"{raw_tick:.2f}"
        ty = by + 4   # a few px below the top edge so it sits nicely
        g.append(txt(bx + box_w + tick_gap, ty, tick_label, fs=font_size - 1))

    # Bottom tick (value at the very bottom = raw_bottom, displayed as its value)
    bottom_y = boxes_y0 + n_boxes * box_h + 4
    bottom_label = f"{raw_bottom:.2f}"
    g.append(txt(x0 + box_w + tick_gap, bottom_y, bottom_label, fs=font_size - 1))

    # ── Append legend group to SVG root ──────────────────────────────────────
    root.append(g)


def add_gene_title(root, gene_name, font_size=18, padding_bottom=10):
    """
    Add the gene name to the SVG as:
      1. A <title> element (SVG metadata — shown as tooltip in browsers)
      2. A visible <text> element centred at the top of the canvas

    The existing SVG content is wrapped in a <g> and shifted down to make
    room for the title text plus padding_bottom spacing beneath it.

    Args:
        root: the root <svg> ElementTree element (modified in-place)
        gene_name: string to display
        font_size: font size in px (default 18)
        padding_bottom: extra space in SVG units between the title baseline
                        and the original content (default 10)
    """
    SVG_NS = 'http://www.w3.org/2000/svg'

    # Total vertical offset = font size + padding below the title
    offset_y = font_size + padding_bottom

    # ── 1. <title> metadata element ──────────────────────────────────────────
    existing_title = root.find(f'{{{SVG_NS}}}title')
    if existing_title is not None:
        root.remove(existing_title)

    title_elem = ET.Element(f'{{{SVG_NS}}}title')
    title_elem.text = gene_name
    root.insert(0, title_elem)

    # ── 2. Expand canvas height to fit the title ─────────────────────────────
    height_attr = root.get('height', '')
    if height_attr:
        try:
            new_h = float(''.join(c for c in height_attr if c.isdigit() or c == '.')) + offset_y
            unit = ''.join(c for c in height_attr if c.isalpha())
            root.set('height', f"{new_h}{unit}")
        except ValueError:
            pass

    vb = root.get('viewBox', '')
    if vb:
        try:
            parts = vb.split()
            parts[3] = str(float(parts[3]) + offset_y)
            root.set('viewBox', ' '.join(parts))
        except (IndexError, ValueError):
            pass

    # ── 3. Wrap all existing content in a <g> shifted down by offset_y ───────
    children = [c for c in list(root) if c is not title_elem]
    for child in children:
        root.remove(child)

    wrapper = ET.Element(f'{{{SVG_NS}}}g')
    wrapper.set('transform', f'translate(0, {offset_y})')
    for child in children:
        wrapper.append(child)
    root.append(wrapper)

    # ── 4. Visible <text> label centred at the top ───────────────────────────
    width_attr = root.get('width', '')
    cx = None
    if width_attr:
        try:
            cx = float(''.join(c for c in width_attr if c.isdigit() or c == '.')) / 2
        except ValueError:
            pass
    if cx is None:
        try:
            cx = float(vb.split()[2]) / 2
        except (IndexError, ValueError):
            cx = 300  # fallback

    text_elem = ET.Element(f'{{{SVG_NS}}}text')
    text_elem.set('x', str(cx))
    text_elem.set('y', str(font_size))   # baseline sits at font_size px from top
    text_elem.set('text-anchor', 'middle')
    text_elem.set('dominant-baseline', 'auto')
    text_elem.set('font-family', 'Arial, sans-serif')
    text_elem.set('font-size', str(font_size))
    text_elem.set('font-weight', 'bold')
    text_elem.set('fill', '#222222')
    text_elem.text = gene_name

    # Insert right after <title>, before the wrapper group
    root.insert(1, text_elem)


def color_guard_cells(svg_file, expression_dict, output_file, gene_name=None, global_max=None, opacity=1.0):
    """
    Color guard cells, mesophyll, and vascular circles in SVG based on
    expression data.

    All cell types are colored on a single shared scale:
      - A global shift is applied to make all values non-negative.
      - Every shifted value is divided by global_max (the max single-cell
        expression from the metadata) to get a ratio in [0, 1].
      - Ratio 1.0 → pure red; ratio 0.0 → yellow (yellow_to_red scheme).
      - Negative raw values map to low-yellow after shifting.

    Args:
        svg_file: Input SVG file path
        expression_dict: Dictionary mapping cell_type to expression values
        output_file: Output SVG file path
        gene_name: Gene identifier to embed as SVG title (optional)
        global_max: Max single-cell expression from JSON metadata. Falls back
                    to the max of all average values if None.
        opacity: Fill opacity (0.0 to 1.0)
    """
    tree = ET.parse(svg_file)
    root = tree.getroot()

    conditions = ['D0', 'W0', 'R15', 'W15']
    SVG_NS = 'http://www.w3.org/2000/svg'

    # ── Global shift + scale setup ────────────────────────────────────────────
    # Step 1: shift every value up so the minimum across all cell types is 0.
    # Step 2: find the highest shifted average among the visualized cell types
    #         and use that as the colour ceiling (pure red).
    # The global_max from metadata is stored for reference only.
    shift = compute_global_shift(expression_dict)
    shifted_dict = {k: v + shift for k, v in expression_dict.items()}

    # Colour ceiling = max shifted average across all visualized cell types
    visualized_keys = (
        [f"{c}_Guard"            for c in conditions] +
        [f"{c}_Mesophyll"        for c in conditions] +
        [f"{c}_Epidermal"        for c in conditions] +
        [f"{c}_Trichome"         for c in conditions] +
        [f"{c}_Vascular"         for c in conditions] +
        [f"{c}_Phloem Parenchyma" for c in conditions] +
        [f"{c}_Phloem companion" for c in conditions]
    )
    vis_shifted_values = [shifted_dict[k] for k in visualized_keys if k in shifted_dict]
    colour_ceiling = max(vis_shifted_values) if vis_shifted_values else 1.0

    if global_max is not None:
        print(f"  Global single-cell max (reference): {global_max:.4f}")
    print(f"  Shift applied:    {shift:.4f}")
    print(f"  Colour ceiling:   {colour_ceiling:.4f}  (max shifted avg across visualized cell types → pure red)")

    def get_color(cell_type_key):
        val = shifted_dict.get(cell_type_key)
        if val is None:
            return None
        return expression_to_color(val, colour_ceiling, 'yellow_to_red')

    # ── Guard cells ───────────────────────────────────────────────────────────
    guard_keys = [f"{c}_Guard" for c in conditions]
    guard_display = {k: shifted_dict[k] for k in guard_keys if k in shifted_dict}
    print_color_table("Guard Cell Expression (yellow → red)", guard_display,
                      colour_ceiling, shift, 'yellow_to_red')

    guard_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Guard")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"{condition}_guard":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    guard_colored += 1

    # ── Mesophyll ─────────────────────────────────────────────────────────────
    meso_keys = [f"{c}_Mesophyll" for c in conditions]
    meso_display = {k: shifted_dict[k] for k in meso_keys if k in shifted_dict}
    print_color_table("Mesophyll Expression (yellow → red)", meso_display,
                      colour_ceiling, shift, 'yellow_to_red')

    meso_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Mesophyll")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"mesophyll_{condition}":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    meso_colored += 1

    # ── Epidermal ─────────────────────────────────────────────────────────────
    epid_keys = [f"{c}_Epidermal" for c in conditions]
    epid_display = {k: shifted_dict[k] for k in epid_keys if k in shifted_dict}
    print_color_table("Epidermal Expression (yellow → red)", epid_display,
                      colour_ceiling, shift, 'yellow_to_red')

    epid_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Epidermal")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"{condition}_epidermal":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    epid_colored += 1

    # ── Trichome ──────────────────────────────────────────────────────────────
    trich_keys = [f"{c}_Trichome" for c in conditions]
    trich_display = {k: shifted_dict[k] for k in trich_keys if k in shifted_dict}
    print_color_table("Trichome Expression (yellow → red)", trich_display,
                      colour_ceiling, shift, 'yellow_to_red')

    trich_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Trichome")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"{condition}_trichome":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    trich_colored += 1

    # ── Vascular circle ───────────────────────────────────────────────────────
    vasc_keys = [f"{c}_Vascular" for c in conditions]
    vasc_display = {k: shifted_dict[k] for k in vasc_keys if k in shifted_dict}
    print_color_table("Vascular Expression (yellow → red)", vasc_display,
                      colour_ceiling, shift, 'yellow_to_red')

    vasc_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Vascular")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"vascular_{condition}":
                for circle in elem.iter(f'{{{SVG_NS}}}circle'):
                    circle.set('fill', color)
                    circle.set('fill-opacity', str(opacity))
                    vasc_colored += 1

    # ── Phloem Parenchyma circle ──────────────────────────────────────────────
    pp_keys = [f"{c}_Phloem Parenchyma" for c in conditions]
    pp_display = {k: shifted_dict[k] for k in pp_keys if k in shifted_dict}
    print_color_table("Phloem Parenchyma Expression (yellow → red)", pp_display,
                      colour_ceiling, shift, 'yellow_to_red')

    pp_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Phloem Parenchyma")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"phloem_parenchyma_{condition}" and elem.tag == f'{{{SVG_NS}}}circle':
                elem.set('fill', color)
                elem.set('fill-opacity', str(opacity))
                pp_colored += 1

    # ── Phloem Companion circle ───────────────────────────────────────────────
    pc_keys = [f"{c}_Phloem companion" for c in conditions]
    pc_display = {k: shifted_dict[k] for k in pc_keys if k in shifted_dict}
    print_color_table("Phloem Companion Expression (yellow → red)", pc_display,
                      colour_ceiling, shift, 'yellow_to_red')

    pc_colored = 0
    for condition in conditions:
        color = get_color(f"{condition}_Phloem companion")
        if color is None:
            continue
        for elem in root.iter():
            if elem.get('id') == f"phloem_companion_{condition}" and elem.tag == f'{{{SVG_NS}}}circle':
                elem.set('fill', color)
                elem.set('fill-opacity', str(opacity))
                pc_colored += 1

    # ── Phloem paths: average of Phloem Parenchyma + Phloem companion ─────────
    # Average the two raw values per condition, shift, then colour with colour_ceiling.
    phloem_avg_shifted = {}
    for condition in conditions:
        pp_key = f"{condition}_Phloem Parenchyma"
        pc_key = f"{condition}_Phloem companion"
        if pp_key in expression_dict and pc_key in expression_dict:
            avg_raw = (expression_dict[pp_key] + expression_dict[pc_key]) / 2
            phloem_avg_shifted[condition] = avg_raw + shift

    if phloem_avg_shifted:
        print("=" * 70)
        print("Phloem paths: avg(Phloem Parenchyma, Phloem companion) → yellow to red")
        print("=" * 70)
        print(f"Colour ceiling (shifted): {colour_ceiling:.4f}")
        print("-" * 70)
        print(f"{'Condition':<10} {'PP raw':>10} {'PC raw':>10} {'Avg raw':>10} {'Shifted':>10} {'Ratio':>7} {'Color'}")
        print("-" * 70)
        for condition in conditions:
            if condition in phloem_avg_shifted:
                pp_raw = expression_dict.get(f"{condition}_Phloem Parenchyma", float('nan'))
                pc_raw = expression_dict.get(f"{condition}_Phloem companion", float('nan'))
                avg_raw = (pp_raw + pc_raw) / 2
                sv = phloem_avg_shifted[condition]
                ratio = min(sv / colour_ceiling, 1.0) if colour_ceiling > 0 else 0
                color = expression_to_color(sv, colour_ceiling, 'yellow_to_red')
                print(f"{condition:<10} {pp_raw:>10.4f} {pc_raw:>10.4f} {avg_raw:>10.4f} {sv:>10.4f} {ratio:>7.3f}    {color}")
        print("=" * 70)
        print()

    phloem_colored = 0
    for condition in conditions:
        if condition not in phloem_avg_shifted:
            continue
        color = expression_to_color(phloem_avg_shifted[condition], colour_ceiling, 'yellow_to_red')
        for elem in root.iter():
            if elem.get('id') == f"{condition}_phloem":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    phloem_colored += 1

    # ── Gene name title ───────────────────────────────────────────────────────
    if gene_name:
        add_gene_title(root, gene_name)
        print(f"✓ Gene title added:                {gene_name}")

    # ── Colour scale legend ───────────────────────────────────────────────────
    add_legend(root, gene_name or "Expression", colour_ceiling, shift)
    print(f"✓ Colour legend added (ceiling={colour_ceiling:.4f}, shift={shift:.4f})")

    # ── Write output ──────────────────────────────────────────────────────────
    ET.register_namespace('', SVG_NS)
    tree.write(output_file, encoding='utf-8', xml_declaration=True)
    print(f"✓ Guard cell paths colored:        {guard_colored}")
    print(f"✓ Mesophyll paths colored:         {meso_colored}")
    print(f"✓ Epidermal paths colored:         {epid_colored}")
    print(f"✓ Trichome paths colored:          {trich_colored}")
    print(f"✓ Vascular circles colored:        {vasc_colored}")
    print(f"✓ Phloem parenchyma circles:       {pp_colored}")
    print(f"✓ Phloem companion circles:        {pc_colored}")
    print(f"✓ Phloem paths (averaged) colored: {phloem_colored}")
    print(f"✓ Output written to: {output_file}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Color guard cells, mesophyll, and vascular circles in SVG based on expression data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg
  python3 color_guard_cells_cli.py merged_grid.svg avg_expression.json output.svg --opacity 0.9

Color schemes:
  Guard cells:       yellow → red (scaled within the 4 guard cell values)
  Mesophyll:         yellow → red (scaled within the 4 mesophyll values)
  Vascular (circle): yellow → red (scaled within the 4 vascular values)

Each cell type group is scaled independently so the highest expresser
in each group always gets pure red.
        """
    )

    parser.add_argument('svg_file', help='Input SVG file path')
    parser.add_argument('json_file', help='Input JSON file with expression data')
    parser.add_argument('output_file', help='Output SVG file path')
    parser.add_argument('--opacity', type=float, default=1.0,
                        help='Fill opacity from 0.0 to 1.0 (default: 1.0)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress output messages')

    args = parser.parse_args()

    if not 0.0 <= args.opacity <= 1.0:
        print("Error: Opacity must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)

    try:
        if not args.quiet:
            print("Loading expression data...")
        expression_dict, gene_name, global_max = load_expression_data(args.json_file)

        if not args.quiet:
            print(f"Gene: {gene_name}")
            if global_max is not None:
                print(f"Global max (from metadata): {global_max:.4f}")
            print("Coloring cells...\n")
        color_guard_cells(args.svg_file, expression_dict, args.output_file,
                          gene_name=gene_name, global_max=global_max, opacity=args.opacity)

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
