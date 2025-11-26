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
    # Load h5ad
    adata = ad.read_h5ad(input_h5ad)

    # Check column exists
    if cell_type_column not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_type_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Decide which genes to use
    if gene_list is None:
        gene_list = list(adata.var_names)

    # Convert matrix to dense if sparse
    X = adata.X
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    # Build a DataFrame for convenience
    df = pd.DataFrame(X, columns=adata.var_names, index=adata.obs_names)
    df[cell_type_column] = adata.obs[cell_type_column].values

    # Filter genes
    genes = [g for g in gene_list if g in df.columns]

    # Group by cell type and compute average expression
    avg_df = df.groupby(cell_type_column)[genes].mean()

    # Convert to JSON structure
    output = []
    for cell_type in avg_df.index:
        for gene in genes:
            output.append({
                "cell_type": str(cell_type),
                "gene": gene,
                "average_expression": float(avg_df.loc[cell_type, gene])
            })

    # Write JSON
    with open(output_json, "w") as f:
        json.dump(output, f, indent=4)

    print(f"Saved averaged expression JSON to {output_json}")


# Example usage
if __name__ == "__main__":
    h5ad_to_json_avg(
        input_h5ad="at.h5ad",
        output_json="avg_expression.json",
        cell_type_column="label_majorXcondition",
        gene_list=["AT3G05727"]  # or None for all genes
    )
