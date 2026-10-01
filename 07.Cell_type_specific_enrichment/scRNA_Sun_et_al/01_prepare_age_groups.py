# --------------------- extract Month and expr matrix ---------------------
import scanpy as sc
import pandas as pd
import os
import numpy as np
import scipy.sparse as sp

os.chdir("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT")
ad = sc.read_h5ad("Sun_et_al.h5ad")

#- 1M - 2M
ad_1M_2M = ad[ad.obs['age'].isin(['1M','2M'])].copy()
ad_1M_2M.var_names = ad_1M_2M.var['gene_symbol']
ad_1M_2M.write("Sun/age_1M_2M.h5ad")

#- 5M
ad_5M = ad[ad.obs['age'] == '5M'].copy()
ad_5M.var_names = ad_5M.var['gene_symbol']
ad_5M.write("Sun/age_5M.h5ad")

#- 12M - 15M
ad_12M_15M = ad[ad.obs['age'].isin(['12M','15M'])].copy()
ad_12M_15M.var_names = ad_12M_15M.var['gene_symbol']
ad_12M_15M.write("Sun/age_12M_15M.h5ad")


def mat(adata):
    expr_df = adata.to_df()
    expr_df['cell_type'] = adata.obs['CellType'].values
    
    # per gene avrg expression in per cell type, rm expr 0 in all cell type
    mean_expr_in_each_cell_type = expr_df.groupby('cell_type').mean().T
    mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]

    return(mean_expr_in_each_cell_type)


expr_1M_2M = mat(ad_1M_2M)
expr_1M_2M.to_csv('Sun/expr_1M_2M.txt', sep="\t")

expr_5M = mat(ad_5M)
expr_5M.to_csv('Sun/expr_5M.txt', sep="\t")

expr_12M_15M = mat(ad_12M_15M)
expr_12M_15M.to_csv('Sun/expr_12M_15M.txt', sep="\t")

#================================================================
