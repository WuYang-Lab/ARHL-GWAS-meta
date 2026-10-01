### scDRS 
import scanpy as sc

adata = sc.read_10x_mtx(
  "path/to/10x_dir",     # mtx/tsv目录
  var_names="gene_symbols",
  cache=True
)

adata.var_names_make_unique()
adata.raw = adata.copy()


import scanpy as sc
import pandas as pd 
import re
import os 

os.chdir("f:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT")
df = pd.read_csv("Ranum_et_al_GSE114157_Counts_Matrix.csv", index_col=0)

# cell x genes
adata = sc.AnnData(df.T)
adata.var_names_make_unique()
adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

adata.layers['log1p'] = adata.X.copy()

adata.obs['cell_type'] = pd.Index(adata.obs_names).str.split('_').str[0]
adata.write("Ranum/Ranum.h5ad")




