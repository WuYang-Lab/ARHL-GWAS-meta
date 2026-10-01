### scDRS
import pandas as pd 
import scanpy as sc 

adata = sc.read_h5ad("Jean_et_al.h5ad")
adata.var_names = adata.var['gene_symbol'].astype(str)
adata.var_names_make_unique()

adata.var.drop(columns=['gene_symbol'], inplace=True)
adata = adata[~adata.obs[group_col].isna()].copy()

adata.write("Jean_et_al_SYMBL.h5ad")
