#============================================================#
#
#  scRNA expression treatment of Philippe Jean
#
#============================================================#
# # # 1) generate matrix in python
import pandas as pd 
import scanpy as sc 

SC = sc.read_h5ad("Jean_et_al.h5ad")
expr_df = SC.to_df()
expr_df['cell_type'] = SC.obs['cell_type'].values

# per gene avrg expression in per cell type, rm expr 0 in all cell type
mean_expr_in_each_cell_type = expr_df.groupby('cell_type').mean().T
mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]
mean_expr_in_each_cell_type.to_csv('Sc_P8_12_20.txt', sep="\t")

#================================================================
