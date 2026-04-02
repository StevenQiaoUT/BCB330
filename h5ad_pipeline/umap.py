#!/Library/Frameworks/Python.framework/Versions/3.12/bin/python3
"""
CGI script: Interactive UMAP Visualization

Accepts a gene name and optional h5ad file path via query string, then
streams back a self-contained HTML page with the interactive UMAP plot.

Query parameters:
    gene      : Gene identifier (required), e.g. AT3G05727
    h5ad      : Path to .h5ad file (optional, falls back to H5AD_FILE env var
                or the DEFAULT_H5AD constant below)
    cell_col  : obs column for cell type labels (default: label_v2)

Usage (Apache/lighttpd CGI):
    https://yourserver/cgi-bin/umap_cgi.py?gene=AT3G05727
    https://yourserver/cgi-bin/umap_cgi.py?gene=AT3G05727&h5ad=/data/at.h5ad

Deployment checklist:
    1. Place this file in your server's cgi-bin directory.
    2. chmod +x umap_cgi.py
    3. Ensure the shebang points to the correct Python interpreter.
    4. Set DEFAULT_H5AD below (or the H5AD_FILE environment variable).
    5. The web-server user needs read access to the .h5ad file.
"""

import cgi
import cgitb
import os
import sys
import json
import html
import traceback

# ── Configuration ─────────────────────────────────────────────────────────────
# Hard-code your default .h5ad path here so the script works without
# query-string overrides on your deployment.
DEFAULT_H5AD   = os.environ.get("H5AD_FILE", "/data/at.h5ad")
DEFAULT_COLUMN = "label_major"	#label_v2

# Enable detailed CGI tracebacks in the browser during development.
# Set to False (or remove) in production.
cgitb.enable()


def send_error(message: str, status: str = "400 Bad Request") -> None:
    """Emit a minimal error page and exit."""
    print(f"Status: {status}")
    print("Content-Type: text/html; charset=utf-8")
    print()
    print(f"""<!DOCTYPE html>
<html><head><title>Error</title></head>
<body>
  <h2>Error</h2>
  <pre>{html.escape(message)}</pre>
</body></html>""")
    sys.exit(0)


def build_figure(h5ad_file: str, gene_name: str, cell_col: str):
    """
    Load the h5ad file and build a Plotly Figure object.
    Returns (fig, stats_dict) or raises on failure.
    """
    import scanpy as sc
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np

    adata = sc.read_h5ad(h5ad_file)

    if gene_name not in adata.var_names:
        raise ValueError(
            f"Gene '{gene_name}' not found. "
            f"First 10 available: {list(adata.var_names[:10])}"
        )
    if "X_umap" not in adata.obsm:
        raise ValueError(
            f"No UMAP embedding found. Available: {list(adata.obsm.keys())}"
        )
    if cell_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_col}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Build data frame
    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
    umap_df["cell_type"]       = adata.obs[cell_col].astype(str).values
    umap_df["gene_expression"] = adata[:, gene_name].X.toarray().flatten()

    cell_types = sorted(umap_df["cell_type"].unique())

    # Per-cell-type statistics
    cell_type_stats = {}
    for ct in cell_types:
        expr = umap_df[umap_df["cell_type"] == ct]["gene_expression"]
        cell_type_stats[ct] = {
            "mean":  float(expr.mean()),
            "std":   float(expr.std()),
            "count": int(len(expr)),
        }

    all_mean = float(umap_df["gene_expression"].mean())
    all_std  = float(umap_df["gene_expression"].std())

    # ── Build Plotly figure ───────────────────────────────────────────────────
    fig = go.Figure()

    # Base trace: all cells coloured by expression
    fig.add_trace(go.Scattergl(
        x=umap_df["UMAP_1"],
        y=umap_df["UMAP_2"],
        mode="markers",
        marker=dict(
            size=4,
            color=umap_df["gene_expression"],
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title=f"{gene_name} Expression"),
            opacity=0.7,
        ),
        customdata=np.column_stack([umap_df["cell_type"], umap_df["gene_expression"]]),
        hovertemplate=(
            "<b>Cell Type:</b> %{customdata[0]}<br>"
            "<b>UMAP_1:</b> %{x:.2f}<br>"
            "<b>UMAP_2:</b> %{y:.2f}<br>"
            "<b>Expression:</b> %{customdata[1]:.3f}<br>"
            "<extra></extra>"
        ),
        name="All Cells",
        visible=True,
    ))

    # One highlight trace per cell type
    for ct in cell_types:
        ct_df = umap_df[umap_df["cell_type"] == ct]
        fig.add_trace(go.Scattergl(
            x=ct_df["UMAP_1"],
            y=ct_df["UMAP_2"],
            mode="markers",
            marker=dict(
                size=7,
                color="#FF0000",
                line=dict(width=1, color="black"),
                opacity=1.0,
            ),
            customdata=np.column_stack([ct_df["cell_type"], ct_df["gene_expression"]]),
            hovertemplate=(
                "<b>Cell Type:</b> %{customdata[0]}<br>"
                "<b>UMAP_1:</b> %{x:.2f}<br>"
                "<b>UMAP_2:</b> %{y:.2f}<br>"
                "<b>Expression:</b> %{customdata[1]:.3f}<br>"
                "<extra></extra>"
            ),
            name=ct,
            visible=False,
        ))

    # ── Dropdown buttons ──────────────────────────────────────────────────────
    buttons = [
        dict(
            label=f"All Cells (μ={all_mean:.3f}, σ={all_std:.3f})",
            method="update",
            args=[
                {"visible": [True] + [False] * len(cell_types)},
                {"title": {"text": f"Interactive UMAP - All Cells - {gene_name}"}},
            ],
        )
    ]
    for i, ct in enumerate(cell_types):
        visible = [True] + [False] * len(cell_types)
        visible[i + 1] = True
        stats = cell_type_stats[ct]
        buttons.append(dict(
            label=f"{ct} ({stats['count']:,} cells, μ={stats['mean']:.3f}, σ={stats['std']:.3f})",
            method="update",
            args=[
                {"visible": visible},
                {"title": {"text": f"Interactive UMAP - Highlighting: {ct} - {gene_name}"}},
            ],
        ))

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
            borderwidth=1,
        )],
        title={"text": f"Interactive UMAP - Select Cell Type to Highlight - {gene_name}"},
        xaxis_title="UMAP_1",
        yaxis_title="UMAP_2",
        template="plotly_white",
        showlegend=False,
    )

    stats_out = {
        "gene":       gene_name,
        "total_cells": int(len(umap_df)),
        "cell_types": len(cell_types),
        "all_mean":   all_mean,
        "all_std":    all_std,
        "per_cell_type": cell_type_stats,
    }
    return fig, stats_out


def main() -> None:
    # ── Parse query string ────────────────────────────────────────────────────
    form     = cgi.FieldStorage()
    gene     = form.getfirst("gene", "").strip()
    h5ad     = form.getfirst("h5ad",  DEFAULT_H5AD).strip()
    cell_col = form.getfirst("cell_col", DEFAULT_COLUMN).strip()

    # ── Validate ──────────────────────────────────────────────────────────────
    if not gene:
        send_error(
            "Missing required parameter: gene\n\n"
            "Usage: umap_cgi.py?gene=AT3G05727"
        )

    if not os.path.isfile(h5ad):
        send_error(
            f"h5ad file not found: {h5ad}\n"
            "Set DEFAULT_H5AD in the script or pass ?h5ad=/path/to/file.h5ad",
            status="500 Internal Server Error",
        )

    # ── Build plot ────────────────────────────────────────────────────────────
    try:
        fig, stats = build_figure(h5ad, gene, cell_col)
    except ValueError as exc:
        send_error(str(exc))
    except Exception:
        send_error(
            "Unexpected error while building plot:\n\n" + traceback.format_exc(),
            status="500 Internal Server Error",
        )

    # ── Stream self-contained HTML ────────────────────────────────────────────
    # fig.to_html() returns a full page; we add a small stats footer.
    html_body = fig.to_html(
        full_html=True,
        include_plotlyjs="cdn",   # keeps response small; use "inline" if offline
        config={"responsive": True},
    )

    # Inject a lightweight stats block before </body>
    stats_html = f"""
<div style="font-family:monospace;font-size:12px;padding:8px 16px;
            background:#f5f5f5;border-top:1px solid #ddd;">
  <b>Gene:</b> {html.escape(stats['gene'])} &nbsp;|&nbsp;
  <b>Total cells:</b> {stats['total_cells']:,} &nbsp;|&nbsp;
  <b>Cell types:</b> {stats['cell_types']} &nbsp;|&nbsp;
  <b>Mean expr (all):</b> {stats['all_mean']:.4f} &nbsp;|&nbsp;
  <b>Std:</b> {stats['all_std']:.4f}
</div>
"""
    html_body = html_body.replace("</body>", stats_html + "\n</body>", 1)

    print("Status: 200 OK")
    print("Content-Type: text/html; charset=utf-8")
    print()
    print(html_body)


if __name__ == "__main__":
    main()
