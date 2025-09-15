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
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

dir <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Iyer"
dt <- fread(file.path(dir, "expr.txt"))

dt <- dt %>% gather(celltype, exp_norm, -gene)
dt$celltype <- gsub("[[:space:]/-]", "_", dt$celltype)
dt$celltype <- gsub("'", "_", dt$celltype)

#- keep only mouse to human 1to1 mapped orthologs
conf <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision"
m2h <- fread(file.path(conf, "mouse_human_homologs.txt"))
names(m2h) <- c("MumSYM", "HumSYM")  

dat <- merge(dt, m2h, by.x="gene", by.y="MumSYM")


#- fill in NAs with 0:
tmp <- dat %>% pivot_wider(id_cols=c(gene, HumSYM), names_from=celltype, values_from=exp_norm, values_fill=0, values_fn=list(exp_norm=sum))
dat <- tmp %>% pivot_longer(cols=-c(gene, HumSYM), names_to="celltype", values_to="exp_norm")

#- remove duplicated genes
tmp.dup <- dat %>% add_count(gene) 

#- remove genes not expressed in any cell type
tmp.noexp <- dat %>% 
  group_by(gene) %>% 
  summarise(sum_exp=sum(exp_norm)) %>%
  filter(sum_exp==0)

dat <- dat %>% filter(!gene %in% tmp.noexp$gene)

#- Calculate specificity  
dat <- dat %>%
  group_by(gene) %>%
  mutate(specificity=exp_norm/sum(exp_norm)) %>%
  ungroup()

## Keep only genes tested in MAGMA  
# 4. Write MAGMA and LDSC input  
#  Get number of genes that represent 10% of the dataset
n_genes <- length(unique(dat$HumSYM))
n_genes_to_keep <- (n_genes * 0.1) %>% round()

### Get MAGMA input top10%
magma_top10 <- function(d, Cell_type){
  d_spe <- d %>% group_by(.data[[Cell_type]]) %>% slice_max(order_by=specificity, n=n_genes_to_keep, with_ties=TRUE)
  d_spe %>% group_split(.keep = TRUE) %>% purrr::walk(~ write_group_magma(.x, Cell_type))
}

write_group_magma <- function(df,Cell_type) {
  df <- dplyr::select(df, all_of(Cell_type), HumSYM)
  df_name <- make.names(unique(df[Cell_type]))
  colnames(df)[2] <- df_name  
  dir.create("MAGMA", showWarnings = FALSE)
  
  dplyr::select(df, 2) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Cat") %>%
    readr::write_tsv("MAGMA/top10.txt", append=TRUE)
  invisible(df)
}

### Get LDSC input top 10%
write_group  = function(df,Cell_type) {
  df <- dplyr::select(df, dplyr::all_of(Cell_type), HumSYM)
  
  dir.create("LDSC", showWarnings=FALSE)   
  df %>% dplyr::select(HumSYM) %>% readr::write_tsv(paste0("LDSC/",make.names(unique(df[1])),".bed"), col_names=F)
  invisible(df)
}

ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by(.data[[Cell_type]]) %>% slice_max(order_by=specificity, n=n_genes_to_keep, with_ties=TRUE)
  d_spe %>% group_split(.keep = TRUE) %>% purrr::walk(~ write_group(.x, Cell_type))
}

### Write MAGMA/LDSC input files 
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Iyer")
if (file.exists("MAGMA/top10.txt")) {file.remove("MAGMA/top10.txt")}
dat %>% filter(exp_norm>1) %>% magma_top10("celltype")

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Iyer")
dat %>% filter(exp_norm>1) %>% ldsc_bedfile("celltype")
control <- as.data.frame(unique(dat$HumSYM))
readr::write_tsv(control,  "LDSC/control.bed", col_names=F)