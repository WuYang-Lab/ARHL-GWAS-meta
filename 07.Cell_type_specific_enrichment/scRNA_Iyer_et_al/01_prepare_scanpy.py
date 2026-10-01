#============================================================#
#
#  scRNA expression treatment of Iyer et al
#
#============================================================#
import scanpy as sc
import pandas as pd
import os
import numpy as np


os.chdir("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT")
ad = sc.read_h5ad("Iyer_et_al.h5ad")

ad.obs['celltype'] = ad.obs['clusters']
sub_ad = ad[ad.obs['celltype'] != 'unknown'].copy()

sub_ad.var_names = sub_ad.var['gene_symbol']
sub_ad.write("Iyer/Iyer_et_al.h5ad")

expr_df = sub_ad.to_df()
expr_df['cell_type'] = adata.obs['CellType'].values
# per gene avrg expression in per cell type, rm expr 0 in all cell type
mean_expr_in_each_cell_type = expr_df.groupby('cell_type').mean().T
mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]

mean_expr_in_each_cell_type.to_csv('Iyer/expr.txt', sep="\t")

#================================================================
