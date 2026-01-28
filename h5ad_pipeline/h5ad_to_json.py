"""
Convert H5AD to JSON with Average Expression by Cell Type

Computes average gene expression for each cell type and exports to JSON format.

Usage:
    python h5ad_to_json_avg.py <input_h5ad> <output_json> <cell_type_column> [gene1 gene2 ...]

Arguments:
    input_h5ad        : Path to input .h5ad file
    output_json       : Path to output JSON file
    cell_type_column  : Column name in adata.obs containing cell type labels
    gene1 gene2 ...   : (Optional) List of genes to include. If not provided, uses all genes.

Examples:
    # Single gene
    python h5ad_to_json_avg.py at.h5ad output.json label_majorXcondition AT3G05727

    # Multiple genes
    python h5ad_to_json_avg.py at.h5ad output.json label_majorXcondition AT3G05727 AT1G01010 AT5G12345

    # All genes (no gene list specified)
    python h5ad_to_json_avg.py at.h5ad output.json label_v2

    # Different paths
    python h5ad_to_json_avg.py /path/to/data.h5ad /path/to/output.json cell_type
"""

import sys
import os
import json
import anndata as ad
import numpy as np
import pandas as pd


def h5ad_to_json_avg(
        input_h5ad,
        output_json,
        cell_type_column="cell_type",
        gene_list=None
):
    """
    Convert h5ad file to JSON with averaged expression per cell type.

    Parameters
    ----------
    input_h5ad : str
        Path to input .h5ad file
    output_json : str
        Path to output JSON file
    cell_type_column : str
        Column name in adata.obs containing cell type labels
    gene_list : list of str or None
        List of genes to include. If None, uses all genes.
    """
    # Load h5ad
    print(f"Loading data from: {input_h5ad}")
    try:
        adata = ad.read_h5ad(input_h5ad)
    except Exception as e:
        print(f"Error reading h5ad file: {e}")
        sys.exit(1)

    print(f"  Total cells: {adata.n_obs:,}")
    print(f"  Total genes: {adata.n_vars:,}")

    # Check column exists
    if cell_type_column not in adata.obs.columns:
        print(f"\nError: Column '{cell_type_column}' not found in dataset")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)

    # Get unique cell types
    unique_cell_types = adata.obs[cell_type_column].unique()
    print(f"  Cell types in '{cell_type_column}': {len(unique_cell_types)}")

    # Decide which genes to use
    if gene_list is None:
        gene_list = list(adata.var_names)
        print(f"  Using all {len(gene_list):,} genes")
    else:
        # Validate genes exist
        missing_genes = [g for g in gene_list if g not in adata.var_names]
        if missing_genes:
            print(f"\nWarning: The following genes were not found in dataset:")
            for g in missing_genes:
                print(f"  - {g}")
            print(f"\nAvailable genes (first 10): {list(adata.var_names[:10])}")

        # Filter to only existing genes
        gene_list = [g for g in gene_list if g in adata.var_names]
        if len(gene_list) == 0:
            print("\nError: No valid genes specified")
            sys.exit(1)
        print(f"  Using {len(gene_list)} gene(s): {gene_list}")

    # Convert matrix to dense if sparse
    print("\nProcessing expression data...")
    X = adata.X
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    # Build a DataFrame for convenience
    df = pd.DataFrame(X, columns=adata.var_names, index=adata.obs_names)
    df[cell_type_column] = adata.obs[cell_type_column].values

    # Filter genes
    genes = [g for g in gene_list if g in df.columns]

    # Group by cell type and compute average expression
    print("Computing average expression per cell type...")
    avg_df = df.groupby(cell_type_column)[genes].mean()

    print(f"  Computed averages for {len(avg_df)} cell types")

    # Convert to JSON structure
    print("Converting to JSON format...")
    output = []
    for cell_type in avg_df.index:
        for gene in genes:
            output.append({
                "cell_type": str(cell_type),
                "gene": gene,
                "average_expression": float(avg_df.loc[cell_type, gene])
            })

    # Write JSON
    try:
        with open(output_json, "w") as f:
            json.dump(output, f, indent=4)
        print(f"\n✓ Successfully saved to: {output_json}")
        print(f"  Total records: {len(output):,}")
    except Exception as e:
        print(f"\nError writing JSON file: {e}")
        sys.exit(1)


def main():
    # Parse command-line arguments
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    input_h5ad = sys.argv[1]
    output_json = sys.argv[2]
    cell_type_column = sys.argv[3]
    gene_list = sys.argv[4:] if len(sys.argv) > 4 else None

    # Validate input file exists
    if not os.path.exists(input_h5ad):
        print(f"Error: Input file '{input_h5ad}' does not exist")
        sys.exit(1)

    if not input_h5ad.endswith('.h5ad'):
        print(f"Warning: File '{input_h5ad}' does not have .h5ad extension")

    # Print configuration
    print("=" * 60)
    print("H5AD to JSON Converter - Average Expression")
    print("=" * 60)
    print(f"Input file:       {input_h5ad}")
    print(f"Output file:      {output_json}")
    print(f"Cell type column: {cell_type_column}")
    if gene_list:
        print(f"Gene list:        {', '.join(gene_list)}")
    else:
        print(f"Gene list:        All genes")
    print("=" * 60)
    print()

    # Run conversion
    h5ad_to_json_avg(
        input_h5ad=input_h5ad,
        output_json=output_json,
        cell_type_column=cell_type_column,
        gene_list=gene_list
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
