#!/mnt/home/sqiao/venv/bin/python3
"""
eFP Viewer — CGI Script
========================
Loads an H5AD file once, then produces a single HTML page containing:
  1. A colored eFP SVG (tissue expression browser)
  2. An interactive Plotly UMAP with cell-type dropdown

URL Parameters (minimal — dataset and column are chosen via in-page dropdowns):
    gene : Gene ID (e.g. AT3G05727)
    ds   : Dataset key — "allcells" or "selectcells" (default: allcells)
    col  : obs column key (default: label_majorXcondition)

All other parameters (file path, SVG template, umapcol, opacity)
are resolved server-side from the dataset key and are never exposed in the URL.

Usage (CGI / browser):
    http://142.150.215.219/~sqiao/cgi-bin/efp_viewer.py?gene=AT3G05727&ds=allcells&col=label_majorXcondition

Usage (command-line):
    python efp_viewer.py <ds_key> <gene> [output.html]
    python efp_viewer.py allcells AT3G05727 out.html
"""

import sys
import os
import io
import warnings
from xml.etree import ElementTree as ET

warnings.filterwarnings("ignore", category=FutureWarning)


# ──────────────────────────────────────────────────────────────────────────────
# Dataset registry — all server-side config lives here, never in the URL
# ──────────────────────────────────────────────────────────────────────────────

DATASETS = {
    "allcells": {
        "label":    "All Cells",
        "h5ad":     "/mnt/home/sqiao/h5ad_files/allcells.h5ad",
        "svg":      "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/svg_template.svg",
        "umap_col": "label_major",
        "opacity":  1.0,
    },
    "selectcells": {
        "label":    "Select Cells",
        "h5ad":     "/mnt/home/sqiao/h5ad_files/selectedcells.h5ad",
        "svg":      "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/svg_template.svg",
        "umap_col": "label_major",
        "opacity":  1.0,
    },
}

# Column registry — maps URL key → display label shown in the dropdown.
# Add new entries here when future h5ad types introduce additional obs columns.
COLUMNS = {
    "label_majorXcondition": "label_majorXcondition",
}

DEFAULT_DATASET = "allcells"
DEFAULT_GENE    = "SZF1"
DEFAULT_COLUMN  = "label_majorXcondition"


# ──────────────────────────────────────────────────────────────────────────────
# Color utilities
# ──────────────────────────────────────────────────────────────────────────────

def expression_to_color(value, max_val, scheme='yellow_to_red'):
    if scheme == 'yellow_to_red':
        minColor = {'red': 255, 'green': 255, 'blue': 0}
    else:
        minColor = {'red': 255, 'green': 255, 'blue': 255}
    maxColor = {'red': 255, 'green': 0, 'blue': 0}
    ratio = min(value / max_val, 1.0) if max_val > 0 else 0
    if not (0 <= ratio <= 1):
        return "#ffffff"
    red   = minColor['red']   + round((maxColor['red']   - minColor['red'])   * ratio)
    green = minColor['green'] + round((maxColor['green'] - minColor['green']) * ratio)
    blue  = minColor['blue']  + round((maxColor['blue']  - minColor['blue'])  * ratio)
    return f"#{red:02x}{green:02x}{blue:02x}"


def compute_global_shift(expression_dict):
    vals = list(expression_dict.values())
    mn = min(vals) if vals else 0
    return -mn if mn < 0 else 0


# ──────────────────────────────────────────────────────────────────────────────
# SVG legend + title
# ──────────────────────────────────────────────────────────────────────────────

def add_legend(root, gene_name, colour_ceiling, shift, n_boxes=10):
    SVG_NS = "http://www.w3.org/2000/svg"
    box_w = 24
    box_h = 20
    tick_gap = 6
    font_size = 10
    header_fs = 11
    left_pad = 10
    tick_label_w = 40
    right_pad = 20
    legend_width = left_pad + box_w + tick_gap + tick_label_w + right_pad

    for attr in ["width", "viewBox"]:
        val = root.get(attr, "")
        if not val:
            continue
        try:
            if attr == "width":
                num = float("".join(c for c in val if c.isdigit() or c == "."))
                unit = "".join(c for c in val if c.isalpha())
                root.set(attr, f"{num + legend_width}{unit}")
            else:
                parts = val.split()
                parts[2] = str(float(parts[2]) + legend_width)
                root.set(attr, " ".join(parts))
        except (ValueError, IndexError):
            pass

    existing = [c for c in list(root) if c.get("id") != "expression_legend"]
    for child in existing:
        root.remove(child)

    shift_g = ET.Element(f"{{{SVG_NS}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    x0 = left_pad
    try:
        h_attr = root.get("height", "")
        canvas_h = float("".join(c for c in h_attr if c.isdigit() or c == "."))
    except ValueError:
        canvas_h = 600

    header_h = 3 * (header_fs + 4)
    legend_h = header_h + n_boxes * box_h + font_size + 6
    y0 = canvas_h - legend_h - 10

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

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

    hy = y0 + header_fs
    g.append(txt(x0, hy, gene_name, fs=header_fs, weight="bold", style="italic"))
    hy += header_fs + 4
    g.append(txt(x0, hy, "Single Cell Max"))
    hy += font_size + 4
    g.append(txt(x0, hy, "Linear"))

    boxes_y0 = y0 + header_h
    raw_top    = colour_ceiling - shift
    raw_bottom = -shift

    for i in range(n_boxes):
        ratio_mid = 1.0 - (i + 0.5) / n_boxes
        color = expression_to_color(ratio_mid * colour_ceiling, colour_ceiling, "yellow_to_red")
        by = boxes_y0 + i * box_h
        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x", str(x0))
        rect.set("y", str(by))
        rect.set("width", str(box_w))
        rect.set("height", str(box_h))
        rect.set("fill", color)
        rect.set("stroke", "#888888")
        rect.set("stroke-width", "0.5")
        g.append(rect)

        ratio_top = 1.0 - i / n_boxes
        raw_tick = raw_bottom + ratio_top * (raw_top - raw_bottom)
        g.append(txt(x0 + box_w + tick_gap, by + 4, f"{raw_tick:.2f}", fs=font_size - 1))

    g.append(txt(x0 + box_w + tick_gap, boxes_y0 + n_boxes * box_h + 4,
                 f"{raw_bottom:.2f}", fs=font_size - 1))
    root.append(g)


def add_gene_title(root, gene_name, font_size=18, padding_bottom=10):
    SVG_NS = 'http://www.w3.org/2000/svg'
    offset_y = font_size + padding_bottom

    existing = root.find(f'{{{SVG_NS}}}title')
    if existing is not None:
        root.remove(existing)

    t = ET.Element(f'{{{SVG_NS}}}title')
    t.text = gene_name
    root.insert(0, t)

    for attr in ['height', 'viewBox']:
        val = root.get(attr, '')
        if not val:
            continue
        try:
            if attr == 'height':
                num  = float(''.join(c for c in val if c.isdigit() or c == '.')) + offset_y
                unit = ''.join(c for c in val if c.isalpha())
                root.set(attr, f"{num}{unit}")
            else:
                parts = val.split()
                parts[3] = str(float(parts[3]) + offset_y)
                root.set(attr, ' '.join(parts))
        except (ValueError, IndexError):
            pass

    children = [c for c in list(root) if c is not t]
    for c in children:
        root.remove(c)

    wrap = ET.Element(f'{{{SVG_NS}}}g')
    wrap.set('transform', f'translate(0, {offset_y})')
    for c in children:
        wrap.append(c)
    root.append(wrap)

    try:
        cx = float(''.join(c for c in root.get('width', '') if c.isdigit() or c == '.')) / 2
    except ValueError:
        cx = 300

    tx = ET.Element(f'{{{SVG_NS}}}text')
    tx.set('x', str(cx))
    tx.set('y', str(font_size))
    tx.set('text-anchor', 'middle')
    tx.set('font-family', 'Arial, sans-serif')
    tx.set('font-size', str(font_size))
    tx.set('font-weight', 'bold')
    tx.set('fill', '#222222')
    tx.text = gene_name
    root.insert(1, tx)


# ──────────────────────────────────────────────────────────────────────────────
# SVG coloring
# ──────────────────────────────────────────────────────────────────────────────

def color_svg(svg_file, expression_dict, gene_name=None, opacity=1.0):
    tree = ET.parse(svg_file)
    root = tree.getroot()
    SVG_NS = 'http://www.w3.org/2000/svg'
    conditions = ['D0', 'W0', 'R15', 'W15']

    shift = compute_global_shift(expression_dict)
    shifted_dict = {k: v + shift for k, v in expression_dict.items()}

    vis_keys = (
        [f"{c}_Guard"             for c in conditions] +
        [f"{c}_Mesophyll"         for c in conditions] +
        [f"{c}_Epidermal"         for c in conditions] +
        [f"{c}_Trichome"          for c in conditions] +
        [f"{c}_Vascular"          for c in conditions] +
        [f"{c}_Phloem Parenchyma" for c in conditions] +
        [f"{c}_Phloem companion"  for c in conditions]
    )
    vis_vals = [shifted_dict[k] for k in vis_keys if k in shifted_dict]
    colour_ceiling = max(vis_vals) if vis_vals else 1.0

    def get_color(key):
        v = shifted_dict.get(key)
        return expression_to_color(v, colour_ceiling, 'yellow_to_red') if v is not None else None

    def raw_val(key):
        v = expression_dict.get(key)
        return f"{v:.4f}" if v is not None else "N/A"

    def add_tooltip(elem, text):
        t = ET.Element(f'{{{SVG_NS}}}title')
        t.text = text
        elem.insert(0, t)

    # Guard
    for c in conditions:
        key = f"{c}_Guard"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"{c}_guard":
                for p in el.iter(f'{{{SVG_NS}}}path'):
                    p.set('fill', color)
                    p.set('fill-opacity', str(opacity))
                    add_tooltip(p, tip)

    # Mesophyll
    for c in conditions:
        key = f"{c}_Mesophyll"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"mesophyll_{c}":
                for p in el.iter(f'{{{SVG_NS}}}path'):
                    p.set('fill', color)
                    p.set('fill-opacity', str(opacity))
                    add_tooltip(p, tip)

    # Epidermal
    for c in conditions:
        key = f"{c}_Epidermal"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"{c}_epidermal":
                for p in el.iter(f'{{{SVG_NS}}}path'):
                    p.set('fill', color)
                    p.set('fill-opacity', str(opacity))
                    add_tooltip(p, tip)

    # Trichome
    for c in conditions:
        key = f"{c}_Trichome"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"{c}_trichome":
                for p in el.iter(f'{{{SVG_NS}}}path'):
                    p.set('fill', color)
                    p.set('fill-opacity', str(opacity))
                    add_tooltip(p, tip)

    # Vascular
    for c in conditions:
        key = f"{c}_Vascular"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"vascular_{c}":
                for circle in el.iter(f'{{{SVG_NS}}}circle'):
                    circle.set('fill', color)
                    circle.set('fill-opacity', str(opacity))
                    add_tooltip(circle, tip)

    # Phloem Parenchyma
    for c in conditions:
        key = f"{c}_Phloem Parenchyma"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"phloem_parenchyma_{c}" and el.tag == f'{{{SVG_NS}}}circle':
                el.set('fill', color)
                el.set('fill-opacity', str(opacity))
                add_tooltip(el, tip)

    # Phloem Companion
    for c in conditions:
        key = f"{c}_Phloem companion"
        color = get_color(key)
        if not color:
            continue
        tip = f"{key}\nAvg expression: {raw_val(key)}"
        for el in root.iter():
            if el.get('id') == f"phloem_companion_{c}" and el.tag == f'{{{SVG_NS}}}circle':
                el.set('fill', color)
                el.set('fill-opacity', str(opacity))
                add_tooltip(el, tip)

    # Phloem paths (averaged)
    for c in conditions:
        pp_key = f"{c}_Phloem Parenchyma"
        pc_key = f"{c}_Phloem companion"
        if pp_key not in expression_dict or pc_key not in expression_dict:
            continue

        pp_raw  = expression_dict[pp_key]
        pc_raw  = expression_dict[pc_key]
        avg_raw = (pp_raw + pc_raw) / 2
        color   = expression_to_color(avg_raw + shift, colour_ceiling, 'yellow_to_red')
        tip = (
            f"{c}_Phloem (averaged)\n"
            f"Phloem Parenchyma: {pp_raw:.4f}\n"
            f"Phloem Companion: {pc_raw:.4f}\n"
            f"Avg expression: {avg_raw:.4f}"
        )
        for el in root.iter():
            if el.get('id') == f"{c}_phloem":
                for p in el.iter(f'{{{SVG_NS}}}path'):
                    p.set('fill', color)
                    p.set('fill-opacity', str(opacity))
                    add_tooltip(p, tip)

    if gene_name:
        add_gene_title(root, gene_name)
    add_legend(root, gene_name or "Expression", colour_ceiling, shift)
    return tree


# ──────────────────────────────────────────────────────────────────────────────
# H5AD loader
# ──────────────────────────────────────────────────────────────────────────────

def load_all(h5ad_path, cell_type_column, gene_list, umap_col='label_major'):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    print(f"Loading {h5ad_path}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes", file=sys.stderr)

    if cell_type_column not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_type_column}' not found. Available: {list(adata.obs.columns)}"
        )

    # ── gene alias resolution ──────────────────────────────────────────────
    # Build two-way lookup: symbol ↔ TAIR ID
    tair_col = None
    for candidate in ("TAIR_ID", "tair_id", "gene_id"):
        if candidate in adata.var.columns:
            tair_col = candidate
            break

    sym_to_tair = {}
    tair_to_sym = {}
    if tair_col:
        for sym, tair in zip(adata.var_names, adata.var[tair_col]):
            if isinstance(tair, str) and tair:
                sym_to_tair[sym.upper()]  = tair.upper()
                tair_to_sym[tair.upper()] = sym
        print(f"  Alias column found: '{tair_col}' ({len(sym_to_tair)} mappings)", file=sys.stderr)
    else:
        print("  No TAIR_ID alias column found — symbol-only lookup", file=sys.stderr)

    # Resolve each input query to an actual var_name
    var_names_upper = {v.upper(): v for v in adata.var_names}

    def resolve(query):
        q = query.strip().upper()
        if q in var_names_upper:
            return var_names_upper[q]
        if q in tair_to_sym:
            return tair_to_sym[q]
        return None

    resolved = []
    for g in gene_list:
        r = resolve(g)
        if r:
            resolved.append(r)
        else:
            print(f"Warning: '{g}' not found by symbol or TAIR ID", file=sys.stderr)

    if not resolved:
        raise ValueError("No valid genes found in dataset.")

    gene = resolved[0]
    print(f"  Gene: {gene}", file=sys.stderr)

    # Build display name: "SYMBOL/TAIR_ID" (or just symbol if no alias / same string)
    tair_id = sym_to_tair.get(gene.upper(), "")
    if tair_id and tair_id.upper() != gene.upper():
        display_name = f"{gene}/{tair_id}"
    else:
        display_name = gene

    gene_idx = adata.var_names.get_loc(gene)
    X_col = adata.X[:, gene_idx]

    if sp.issparse(X_col):
        X_col = X_col.toarray().ravel()
    else:
        X_col = np.asarray(X_col).ravel()

    X_col = X_col.astype(float)

    cell_types        = adata.obs[cell_type_column].astype(str).to_numpy()
    unique_cell_types = pd.unique(cell_types)

    print("Computing per-cell-type averages...", file=sys.stderr)
    expression_dict = {}
    for ct in unique_cell_types:
        mask = (cell_types == ct)
        if mask.sum() == 0:
            continue
        expression_dict[str(ct)] = float(np.mean(X_col[mask]))

    umap_data = None
    if 'X_umap' in adata.obsm:
        umap_col_actual = umap_col if umap_col in adata.obs.columns else cell_type_column

        umap = adata.obsm['X_umap']
        if sp.issparse(umap):
            umap = umap.toarray()
        else:
            umap = np.asarray(umap)

        if umap.ndim != 2 or umap.shape[1] < 2:
            raise ValueError(f"X_umap has invalid shape: {umap.shape}")

        umap   = umap[:, :2].astype(float)
        labels = adata.obs[umap_col_actual].astype(str).to_numpy()

        umap_df = pd.DataFrame({
            'UMAP_1':          umap[:, 0],
            'UMAP_2':          umap[:, 1],
            'cell_type':       labels,
            'gene_expression': X_col
        })

        before  = len(umap_df)
        umap_df = umap_df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=['UMAP_1', 'UMAP_2', 'gene_expression']
        )
        after = len(umap_df)

        print(f"  UMAP rows before filtering: {before:,}", file=sys.stderr)
        print(f"  UMAP rows after filtering:  {after:,}",  file=sys.stderr)
        print(f"  UMAP coordinates found, using '{umap_col_actual}' for cell labels", file=sys.stderr)

        if after == 0:
            raise ValueError("UMAP data exists, but all rows were filtered out as invalid.")

        umap_data = umap_df
    else:
        print("  Warning: X_umap not found — UMAP panel will be skipped", file=sys.stderr)

    print(f"  Final UMAP dataframe size: {len(umap_data) if umap_data is not None else 0:,}", file=sys.stderr)
    return expression_dict, display_name, umap_data


# ──────────────────────────────────────────────────────────────────────────────
# UMAP Plotly figure → HTML div string
# ──────────────────────────────────────────────────────────────────────────────

def build_umap_div(umap_df, gene_name):
    import plotly.graph_objects as go
    import numpy as np

    if umap_df is None or len(umap_df) == 0:
        return "<div style='padding:20px;color:#b00'>No UMAP points available.</div>"

    x_all      = umap_df['UMAP_1'].to_numpy(dtype=float)
    y_all      = umap_df['UMAP_2'].to_numpy(dtype=float)
    expr_all   = umap_df['gene_expression'].to_numpy(dtype=float)
    labels_all = umap_df['cell_type'].astype(str).to_numpy()

    cell_types = sorted(umap_df['cell_type'].astype(str).unique())
    all_mean   = float(np.mean(expr_all))
    all_std    = float(np.std(expr_all))

    fig = go.Figure()

    fig.add_trace(go.Scattergl(
        x=x_all, y=y_all,
        mode='markers',
        marker=dict(
            size=3,
            color=expr_all,
            colorscale=[[0.0, '#ffff00'], [1.0, '#ff0000']],
            showscale=True,
            colorbar=dict(title=f"{gene_name}<br>Expression"),
            opacity=0.6
        ),
        customdata=np.column_stack([labels_all, expr_all]),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "UMAP1: %{x:.2f}<br>"
            "UMAP2: %{y:.2f}<br>"
            "Expression: %{customdata[1]:.4f}<extra></extra>"
        ),
        name='All Cells',
        visible=True
    ))

    for ct in cell_types:
        ct_mask = (labels_all == ct)
        fig.add_trace(go.Scattergl(
            x=x_all[ct_mask], y=y_all[ct_mask],
            mode='markers',
            marker=dict(
                size=6,
                color='#1a6fcc',
                line=dict(width=0.5, color='#0a3d6b'),
                opacity=1.0
            ),
            customdata=np.column_stack([labels_all[ct_mask], expr_all[ct_mask]]),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "UMAP1: %{x:.2f}<br>"
                "UMAP2: %{y:.2f}<br>"
                "Expression: %{customdata[1]:.4f}<extra></extra>"
            ),
            name=ct,
            visible=False
        ))

    pad   = 0.5
    x_min = float(np.min(x_all)) - pad
    x_max = float(np.max(x_all)) + pad
    y_min = float(np.min(y_all)) - pad
    y_max = float(np.max(y_all)) + pad

    buttons = [{
        "label":  f"All Cells (μ={all_mean:.3f}, σ={all_std:.3f})",
        "method": "update",
        "args": [
            {"visible": [True] + [False] * len(cell_types)},
            {
                "title": {"text": f"UMAP - All Cells - {gene_name}"},
                "xaxis": {"range": [x_min, x_max], "title": "UMAP_1"},
                "yaxis": {"range": [y_min, y_max], "title": "UMAP_2"},
            }
        ]
    }]

    for i, ct in enumerate(cell_types):
        ct_mask  = (labels_all == ct)
        ct_expr  = expr_all[ct_mask]
        ct_n     = int(np.sum(ct_mask))
        ct_mu    = float(np.mean(ct_expr)) if ct_n else float("nan")
        ct_sigma = float(np.std(ct_expr))  if ct_n else float("nan")

        visible = [True] + [False] * len(cell_types)
        visible[i + 1] = True

        buttons.append({
            "label":  f"{ct} (n={ct_n:,}, μ={ct_mu:.3f}, σ={ct_sigma:.3f})",
            "method": "update",
            "args": [
                {"visible": visible},
                {
                    "title": {"text": f"UMAP - {ct} highlighted - {gene_name}"},
                    "xaxis": {"range": [x_min, x_max], "title": "UMAP_1"},
                    "yaxis": {"range": [y_min, y_max], "title": "UMAP_2"},
                }
            ]
        })

    fig.update_layout(
        updatemenus=[dict(
            buttons=buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.01, xanchor="left",
            y=0.99, yanchor="top",
            bgcolor="white",
            bordercolor="gray",
            borderwidth=1
        )],
        title={"text": f"UMAP - All Cells - {gene_name}"},
        xaxis=dict(title="UMAP_1", range=[x_min, x_max]),
        yaxis=dict(title="UMAP_2", range=[y_min, y_max]),
        template="plotly_white",
        showlegend=False,
        height=650,
        margin=dict(t=60, b=40, l=40, r=40)
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn')


# ──────────────────────────────────────────────────────────────────────────────
# HTML page builder
# ──────────────────────────────────────────────────────────────────────────────

def build_html(svg_string, umap_div, gene_name, active_ds_key, active_col):
    """
    active_ds_key : currently-selected dataset key (e.g. "allcells")
    active_col    : currently-selected obs column key (e.g. "label_majorXcondition")

    URL shape produced by this page:
        ?gene=<GENE>&ds=<KEY>&col=<COL>
    No file paths or opacity ever appear in the URL.
    """
    umap_section = ""
    if umap_div:
        umap_section = f"""
        <div class="section-title">Single-Cell UMAP — {gene_name}</div>
        <div id="umap-container">
            {umap_div}
        </div>
"""

    # Build <option> elements for the dataset dropdown
    ds_options = ""
    for key, meta in DATASETS.items():
        sel = ' selected' if key == active_ds_key else ''
        ds_options += f'            <option value="{key}"{sel}>{meta["label"]}</option>\n'

    # Build <option> elements for the column dropdown
    col_options = ""
    for key, label in COLUMNS.items():
        sel = ' selected' if key == active_col else ''
        col_options += f'            <option value="{key}"{sel}>{label}</option>\n'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>eFP Viewer — {gene_name}</title>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: Arial, sans-serif;
            background: #f5f5f5;
            color: #222;
        }}
        header {{
            background: #2c5f2e;
            color: white;
            padding: 16px 32px;
            display: flex;
            align-items: center;
            gap: 24px;
            flex-wrap: wrap;
        }}
        .header-title {{
            font-size: 20px;
            font-weight: bold;
            letter-spacing: 0.5px;
            white-space: nowrap;
        }}
        .header-title span {{
            font-weight: normal;
            font-style: italic;
            margin-left: 12px;
            font-size: 16px;
            opacity: 0.85;
        }}
        .header-controls {{
            display: flex;
            align-items: center;
            gap: 16px;
            margin-left: auto;
            flex-wrap: wrap;
        }}
        .ctrl-group {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .ctrl-group label {{
            font-size: 13px;
            opacity: 0.9;
            white-space: nowrap;
        }}
        /* ── all dropdowns share the same style ── */
        .ctrl-group select {{
            padding: 6px 10px;
            font-size: 14px;
            border: none;
            border-radius: 4px;
            background: white;
            color: #222;
            outline: none;
            cursor: pointer;
            min-width: 130px;
        }}
        .ctrl-group select:focus {{
            box-shadow: 0 0 0 2px #a8d5a2;
        }}
        /* ── gene search ── */
        #gene-input {{
            padding: 6px 10px;
            font-size: 14px;
            border: none;
            border-radius: 4px;
            width: 160px;
            outline: none;
            font-family: monospace;
            letter-spacing: 0.5px;
        }}
        #gene-input:focus {{
            box-shadow: 0 0 0 2px #a8d5a2;
        }}
        #gene-input.invalid {{
            box-shadow: 0 0 0 2px #e05c5c;
            background: #fff5f5;
        }}
        .search-btn {{
            padding: 6px 14px;
            background: #4a8f4c;
            color: white;
            border: none;
            border-radius: 4px;
            font-size: 14px;
            cursor: pointer;
            white-space: nowrap;
            transition: background 0.15s;
        }}
        .search-btn:hover {{ background: #3a7a3c; }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            padding: 24px 32px;
        }}
        .section-title {{
            font-size: 15px;
            font-weight: bold;
            color: #444;
            margin: 28px 0 10px 0;
            padding-bottom: 6px;
            border-bottom: 2px solid #ddd;
        }}
        #svg-container {{
            background: white;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 16px;
            overflow-x: auto;
        }}
        #svg-container svg {{
            display: block;
            max-width: 100%;
            height: auto;
        }}
        #umap-container {{
            background: white;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 8px;
            overflow: hidden;
        }}
    </style>
</head>
<body>
    <header>
        <div class="header-title">
            eFP Browser
            <span>{gene_name}</span>
        </div>

        <form class="header-controls" method="GET" onsubmit="return validateGene()">

            <!-- Dataset picker -->
            <div class="ctrl-group">
                <label for="ds-select">Dataset:</label>
                <select id="ds-select" name="ds">
{ds_options}                </select>
            </div>

            <!-- Column picker -->
            <div class="ctrl-group">
                <label for="col-select">Column:</label>
                <select id="col-select" name="col">
{col_options}                </select>
            </div>

            <!-- Gene search -->
            <div class="ctrl-group">
                <label for="gene-input">Gene ID:</label>
                <input
                    type="text"
                    id="gene-input"
                    name="gene"
                    value="{gene_name.split('/')[0]}"
                    placeholder="e.g. AT3G05727"
                    spellcheck="false"
                    autocomplete="off"
                    autocorrect="off"
                    autocapitalize="characters"
                >
            </div>

            <button type="submit" class="search-btn">&#x1F50D; Search</button>
        </form>
    </header>

    <div class="container">
        <div class="section-title">Tissue eFP — {gene_name}</div>
        <div id="svg-container">
            {svg_string}
        </div>
        {umap_section}
    </div>

</body>
</html>"""


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    is_cgi = "REQUEST_METHOD" in os.environ

    if is_cgi:
        sys.stdout.write("Content-Type: text/html\r\n\r\n")
        sys.stdout.flush()

        from urllib.parse import parse_qs
        params = parse_qs(os.environ.get("QUERY_STRING", ""))

        ds_key    = params.get("ds",   [DEFAULT_DATASET])[0]
        col_key   = params.get("col",  [DEFAULT_COLUMN])[0]
        gene_list = params.get("gene", [DEFAULT_GENE])

        # Unknown keys fall back to defaults
        if ds_key  not in DATASETS: ds_key  = DEFAULT_DATASET
        if col_key not in COLUMNS:  col_key = DEFAULT_COLUMN

        cfg = DATASETS[ds_key]
        h5ad_path        = cfg["h5ad"]
        cell_type_column = col_key          # chosen by the user via dropdown
        svg_template     = cfg["svg"]
        umap_col         = cfg["umap_col"]
        opacity          = cfg["opacity"]

        def html_error(msg):
            sys.stdout.write(
                f"<html><body><pre style='color:red'>Error: {msg}</pre></body></html>"
            )
            sys.stdout.flush()

        if not os.path.exists(h5ad_path):
            html_error(f"H5AD file not found: {h5ad_path}")
            return

        if not os.path.exists(svg_template):
            html_error(f"SVG template not found: {svg_template}")
            return

        try:
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, cell_type_column, gene_list, umap_col
            )

            ET.register_namespace('', 'http://www.w3.org/2000/svg')
            svg_tree = color_svg(
                svg_template,
                expression_dict,
                gene_name=gene_name,
                opacity=opacity
            )

            buf = io.BytesIO()
            svg_tree.write(buf, encoding='utf-8', xml_declaration=False)
            svg_string = buf.getvalue().decode('utf-8')

            umap_div = build_umap_div(umap_df, gene_name) if umap_df is not None else None
            html = build_html(svg_string, umap_div, gene_name,
                              active_ds_key=ds_key, active_col=col_key)

            sys.stdout.write(html)
            sys.stdout.flush()

        except Exception as e:
            import traceback
            print(traceback.format_exc(), file=sys.stderr)
            html_error(str(e))

    else:
        # ── CLI mode ──────────────────────────────────────────────────────────
        import argparse

        parser = argparse.ArgumentParser(description="eFP SVG + UMAP HTML viewer")
        parser.add_argument(
            "ds_key", choices=list(DATASETS.keys()),
            help="Dataset key: " + ", ".join(DATASETS.keys())
        )
        parser.add_argument("gene", nargs='+', help="Gene ID(s); first valid one is used")
        parser.add_argument("output", nargs='?', default="efp_viewer.html")
        parser.add_argument(
            "--col", choices=list(COLUMNS.keys()), default=DEFAULT_COLUMN,
            help="obs column to colour by (default: %(default)s)"
        )
        args = parser.parse_args()

        cfg              = DATASETS[args.ds_key]
        h5ad_path        = cfg["h5ad"]
        cell_type_column = args.col
        svg_template     = cfg["svg"]
        umap_col         = cfg["umap_col"]
        opacity          = cfg["opacity"]

        for f in [h5ad_path, svg_template]:
            if not os.path.exists(f):
                print(f"Error: file not found: {f}")
                sys.exit(1)

        expression_dict, gene_name, umap_df = load_all(
            h5ad_path, cell_type_column, args.gene, umap_col
        )

        ET.register_namespace('', 'http://www.w3.org/2000/svg')
        svg_tree = color_svg(
            svg_template,
            expression_dict,
            gene_name=gene_name,
            opacity=opacity
        )

        buf = io.BytesIO()
        svg_tree.write(buf, encoding='utf-8', xml_declaration=False)
        svg_string = buf.getvalue().decode('utf-8')

        umap_div = build_umap_div(umap_df, gene_name) if umap_df is not None else None
        html = build_html(svg_string, umap_div, gene_name,
                          active_ds_key=args.ds_key, active_col=args.col)

        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(html)

        print(f"✓ Written to: {args.output}")


if __name__ == "__main__":
    main()
