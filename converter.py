import scanpy as sc
import numpy as np
import scipy.sparse as sp
import pandas as pd
from numcodecs import Blosc

# --- load ---
adata = sc.read_h5ad("at.h5ad")

# --- helper: cast to float32 (works for dense & scipy sparse) ---
def to_f32(a):
    if sp.issparse(a):
        return a.astype(np.float32)
    return np.asarray(a, dtype=np.float32)

# --- helper: sanitize dataframe columns conflicting with Zarr writer ---
def sanitize_reserved(df, label=""):
    if "_index" in df.columns:
        # If it's just a duplicate of the actual index, drop it; else rename.
        try:
            same_as_index = df["_index"].astype(str).equals(
                pd.Series(df.index.astype(str), index=df.index)
            )
        except Exception:
            same_as_index = False

        if same_as_index:
            df.drop(columns=["_index"], inplace=True)
            print(f"[sanitize] Dropped redundant '_index' in {label}")
        else:
            new_name = "orig_index" if "orig_index" not in df.columns else "_index_orig"
            df.rename(columns={"_index": new_name}, inplace=True)
            print(f"[sanitize] Renamed '_index' -> '{new_name}' in {label}")

# --- cast main matrix & layers to float32 (optional but often smaller) ---
adata.X = to_f32(adata.X)
for k, v in list(adata.layers.items()):
    adata.layers[k] = to_f32(v)

# cast neighbor graphs (optional)
for k in ("distances", "connectivities"):
    if k in adata.obsp:
        adata.obsp[k] = to_f32(adata.obsp[k])

# --- sanitize obs/var ---
sanitize_reserved(adata.obs, "obs")
sanitize_reserved(adata.var, "var")

# --- sanitize raw.var if raw exists ---
if adata.raw is not None:
    r = adata.raw.to_adata()      # materialize Raw as an AnnData
    sanitize_reserved(r.var, "raw.var")
    # (optional) also make raw.X float32 to save space
    r.X = to_f32(r.X)
    adata.raw = r                  # set back; AnnData will wrap as Raw

# OPTIONAL: make object columns categorical (smaller on disk)
for dfname, df in (("obs", adata.obs), ("var", adata.var)):
    obj_cols = df.select_dtypes(include=["object"]).columns
    for c in obj_cols:
        df[c] = df[c].astype("category")

# --- write with Zstd via Blosc ---
compressor = Blosc(cname="zstd", clevel=7, shuffle=Blosc.BITSHUFFLE)
adata.write_zarr("new_at.zarr")

print("Wrote new_at.zarr successfully.")
