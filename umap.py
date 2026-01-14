import scanpy as sc
import plotly.graph_objects as go
import pandas as pd
import numpy as np

adata = sc.read_h5ad('at.h5ad')

umap_df = pd.DataFrame(
    adata.obsm['X_umap'],
    columns=['UMAP_1', 'UMAP_2']
)

umap_df['cell_type'] = adata.obs['label_v2'].astype(str).values
umap_df['gene_expression'] = adata[:, 'AT3G05727'].X.toarray().flatten()

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
        colorbar=dict(title="Gene Expression"),
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
              {"title": "Interactive UMAP - All Cells"}]
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
                  {"title": f"Interactive UMAP - Highlighting: {cell_type}"}]
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
    title="Interactive UMAP - Select Cell Type to Highlight",
    xaxis_title="UMAP_1",
    yaxis_title="UMAP_2",
    template='plotly_white',
    showlegend=False
)

fig.write_html('interactive_umap_dropdown.html')

print(f"  Total cells: {len(umap_df):,}")
print(f"  Cell types: {len(cell_types)}")
print(f"\nGene Expression Statistics:")
print(f"  All Cells: mean={all_mean:.3f}, std={all_std:.3f}")
for cell_type in cell_types:
    stats = cell_type_stats[cell_type]
    print(f"  {cell_type}: mean={stats['mean']:.3f}, std={stats['std']:.3f}, n={stats['count']:,}")
