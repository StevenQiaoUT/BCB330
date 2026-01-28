"""
Interactive UMAP Visualization with Cell Type Highlighting

Creates an interactive HTML visualization of UMAP embeddings colored by gene expression,
with a dropdown menu to highlight individual cell types.

Usage:
    python interactive_umap.py <h5ad_file> <gene_name> [output_html]

Arguments:
    h5ad_file    : Path to input .h5ad file
    gene_name    : Gene identifier to visualize (e.g., 'AT3G05727')
    output_html  : (Optional) Path to output HTML file (default: interactive_umap.html)

Example:
    python interactive_umap.py at.h5ad AT3G05727 my_umap.html
"""

import sys
import os
import scanpy as sc
import plotly.graph_objects as go
import pandas as pd
import numpy as np


def validate_inputs(h5ad_file, gene_name):
    """Validate input file and gene name"""
    if not os.path.exists(h5ad_file):
        print(f"Error: Input file '{h5ad_file}' does not exist")
        sys.exit(1)

    if not h5ad_file.endswith('.h5ad'):
        print(f"Warning: File '{h5ad_file}' does not have .h5ad extension")


def main():
    # Parse command-line arguments
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    h5ad_file = sys.argv[1]
    gene_name = sys.argv[2]
    output_html = sys.argv[3] if len(sys.argv) > 3 else 'interactive_umap.html'

    # Validate inputs
    validate_inputs(h5ad_file, gene_name)

    print(f"Loading data from: {h5ad_file}")
    print(f"Visualizing gene: {gene_name}")
    print(f"Output will be saved to: {output_html}")
    print()

    # Load data
    try:
        adata = sc.read_h5ad(h5ad_file)
    except Exception as e:
        print(f"Error reading h5ad file: {e}")
        sys.exit(1)

    # Check if gene exists in dataset
    if gene_name not in adata.var_names:
        print(f"Error: Gene '{gene_name}' not found in dataset")
        print(f"Available genes (first 10): {list(adata.var_names[:10])}")
        sys.exit(1)

    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print("Error: UMAP coordinates not found in dataset")
        print("Available embeddings:", list(adata.obsm.keys()))
        sys.exit(1)

    # Check if cell type labels exist
    if 'label_v2' not in adata.obs.columns:
        print("Error: 'label_v2' column not found in dataset")
        print("Available columns:", list(adata.obs.columns))
        sys.exit(1)

    # Extract UMAP coordinates
    umap_df = pd.DataFrame(
        adata.obsm['X_umap'],
        columns=['UMAP_1', 'UMAP_2']
    )

    # Add cell type and gene expression
    umap_df['cell_type'] = adata.obs['label_v2'].astype(str).values
    umap_df['gene_expression'] = adata[:, gene_name].X.toarray().flatten()

    cell_types = sorted(umap_df['cell_type'].unique())

    # Normalize gene expression for color mapping
    expr_min = umap_df['gene_expression'].min()
    expr_max = umap_df['gene_expression'].max()
    umap_df['expr_norm'] = (umap_df['gene_expression'] - expr_min) / (expr_max - expr_min)

    # Calculate statistics for each cell type
    cell_type_stats = {}
    for cell_type in cell_types:
        ct_expr = umap_df[umap_df['cell_type'] == cell_type]['gene_expression']
        cell_type_stats[cell_type] = {
            'mean': ct_expr.mean(),
            'std': ct_expr.std(),
            'count': len(ct_expr)
        }

    # Create figure
    fig = go.Figure()

    # Add base trace (all cells with gene expression color)
    all_mean = umap_df['gene_expression'].mean()
    all_std = umap_df['gene_expression'].std()

    fig.add_trace(go.Scattergl(
        x=umap_df['UMAP_1'],
        y=umap_df['UMAP_2'],
        mode='markers',
        marker=dict(
            size=4,
            color=umap_df['gene_expression'],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title=f"{gene_name} Expression"),
            opacity=0.7
        ),
        customdata=np.column_stack([umap_df['cell_type'], umap_df['gene_expression']]),
        hovertemplate='<b>Cell Type:</b> %{customdata[0]}<br>' +
                      '<b>UMAP_1:</b> %{x:.2f}<br>' +
                      '<b>UMAP_2:</b> %{y:.2f}<br>' +
                      '<b>Gene Expression:</b> %{customdata[1]:.3f}<br>' +
                      '<extra></extra>',
        name='All Cells',
        visible=True
    ))

    # Add a trace for each cell type
    for cell_type in cell_types:
        ct_mask = umap_df['cell_type'] == cell_type
        ct_df = umap_df[ct_mask]

        fig.add_trace(go.Scattergl(
            x=ct_df['UMAP_1'],
            y=ct_df['UMAP_2'],
            mode='markers',
            marker=dict(
                size=7,
                color='#FF0000',
                line=dict(width=1, color='black'),
                opacity=1.0
            ),
            customdata=np.column_stack([ct_df['cell_type'], ct_df['gene_expression']]),
            hovertemplate='<b>Cell Type:</b> %{customdata[0]}<br>' +
                          '<b>UMAP_1:</b> %{x:.2f}<br>' +
                          '<b>UMAP_2:</b> %{y:.2f}<br>' +
                          '<b>Gene Expression:</b> %{customdata[1]:.3f}<br>' +
                          '<extra></extra>',
            name=cell_type,
            visible=False
        ))

    # Create dropdown menu
    buttons = []

    # Button for "All Cells" (show only base trace)
    buttons.append(
        dict(
            label=f"All Cells (μ={all_mean:.3f}, σ={all_std:.3f})",
            method="update",
            args=[{"visible": [True] + [False] * len(cell_types)},
                  {"title": f"Interactive UMAP - All Cells - {gene_name}"}]
        )
    )

    # Button for each cell type
    for i, cell_type in enumerate(cell_types):
        visible = [True] + [False] * len(cell_types)  # Base trace always visible
        visible[i + 1] = True  # Show this cell type trace

        stats = cell_type_stats[cell_type]

        buttons.append(
            dict(
                label=f"{cell_type} ({stats['count']:,} cells, μ={stats['mean']:.3f}, σ={stats['std']:.3f})",
                method="update",
                args=[{"visible": visible},
                      {"title": f"Interactive UMAP - Highlighting: {cell_type} - {gene_name}"}]
            )
        )

    # Update layout with dropdown
    fig.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.01,
                xanchor="left",
                y=0.99,
                yanchor="top",
                bgcolor="white",
                bordercolor="gray",
                borderwidth=1
            )
        ],
        title=f"Interactive UMAP - Select Cell Type to Highlight - {gene_name}",
        xaxis_title="UMAP_1",
        yaxis_title="UMAP_2",
        template='plotly_white',
        showlegend=False
    )

    # Save output
    try:
        fig.write_html(output_html)
        print(f"✓ Successfully saved to: {output_html}")
    except Exception as e:
        print(f"Error saving HTML file: {e}")
        sys.exit(1)

    # Print statistics
    print(f"\n{'=' * 60}")
    print(f"Dataset Statistics")
    print(f"{'=' * 60}")
    print(f"  Input file: {h5ad_file}")
    print(f"  Gene: {gene_name}")
    print(f"  Total cells: {len(umap_df):,}")
    print(f"  Cell types: {len(cell_types)}")
    print(f"\nGene Expression Statistics:")
    print(f"  All Cells: mean={all_mean:.3f}, std={all_std:.3f}")
    print(f"\nPer Cell Type:")
    for cell_type in cell_types:
        stats = cell_type_stats[cell_type]
        print(f"  {cell_type}: mean={stats['mean']:.3f}, std={stats['std']:.3f}, n={stats['count']:,}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
