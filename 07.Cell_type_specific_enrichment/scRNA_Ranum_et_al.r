#============================================================#
#
#  scRNA expression treatment of Ranum
#
#============================================================#
library(dplyr)
library(data.table)
library(stringr)
dir <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT"
data <- fread(file.path(dir, "Ranum_et_al_GSE114157_Counts_Matrix.csv"))

dat_long <- reshape2::melt(data)
colnames(dat_long) <- c("Symbol", "celltype", "exp" )
dat_long <- dat_long %>% separate(celltype, "_", into=c("cell_type_new"), remove=FALSE)
dat <- dat_long %>% group_by(Symbol, cell_type_new) %>% summarise(exp=sum(exp)) 

#- keep only mouse to human 1to1 mapped orthologs
conf <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision"
m2h <- fread(file.path(conf, "mouse_human_homologs.txt"))
names(m2h) <- c("MumSYM", "HumSYM")  

dat <- merge(dat, m2h, by.x="Symbol", by.y="MumSYM")

#- fill in NAs with 0:
tmp <- dat %>% pivot_wider(id_cols=c(Symbol, HumSYM), names_from=cell_type_new, values_from=exp, values_fill=0, values_fn=list(exp=sum))
dat <- tmp %>% pivot_longer(cols=-c(Symbol, HumSYM), names_to="cell_type_new", values_to="exp")

#- remove duplicated genes
tmp.dup <- dat %>% add_count(Symbol) 

#- remove genes not expressed in any cell type
tmp.noexp <- dat %>% group_by(Symbol) %>% summarise(sum_exp=sum(exp)) %>% filter(sum_exp==0)
dat <- dat %>% filter(!Symbol %in% tmp.noexp$Symbol)

# standardize to 1M molecules per cell type
dat <- dat %>% group_by(cell_type_new) %>% mutate(exp_tpm=exp*1e6/sum(exp)) %>% ungroup()

# Calculate specificity  
dat <- dat %>% group_by(Symbol) %>% mutate(specificity=exp_tpm/sum(exp_tpm)) %>% ungroup()
dat <- dat[, c(2:6)]
names(dat) <- c("HumSYM", "CellType", "exp", "exp_tpm", "specificity")

#- Write MAGMA and LDSC input  
#  Get number of genes that represent 10% of the dataset
n_genes <- length(unique(dat$HumSYM))
n_genes_to_keep <- (n_genes * 0.1) %>% round()

### Get MAGMA input top10%
magma_top10 <- function(d,Cell_type){
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

### Get LDSC input top 10% and control
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
# Filter out genes with expression below 1 TPM.
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Ranum/")
if (file.exists("MAGMA/top10.txt")) {file.remove("MAGMA/top10.txt")}
dat %>% filter(exp_tpm>1) %>% magma_top10("CellType")

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Ranum/")
dat %>% filter(exp_tpm>1) %>% ldsc_bedfile("CellType")
control <- as.data.frame(unique(dat$HumSYM))
readr::write_tsv(control,  "LDSC/control.bed", col_names=F)



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




