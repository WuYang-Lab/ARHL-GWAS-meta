#------------------------------------------------------------#
# ---- extract matrix for LDSC-SEG and magma
library(Seurat)
library(tidyverse)
library(data.table)


# Datasets to merge:  
# - stria.rds  
# - stria_im.rds  
# - sgn_rds  
# - sgn_im.rds  

# 1. Normalize cell type per dataset  
##--------------- stria.rds  
Stria <- "/home/wulab/scRNA/scRNA/lu"
# Stria <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_Stria"
stria <- readRDS(file.path(Stria, "stria.rds"))

#- cell ID to cell type
stria.meta <- stria@meta.data %>% select(celltype)
stria.meta$celltype <- gsub(" ","_", stria.meta$celltype)
stria.meta$cell <- rownames(stria.meta)

#- count
stria.counts <- stria@assays$RNA@counts %>% as.data.frame()
stria.counts$gene <- rownames(stria.counts)
stria.counts.long <- stria.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(stria.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(stria.counts.long, file.path(Stria, "stria.counts.long.tsv"), sep="\t", col.names=T)

#- normalized expression
stria.data <- stria@assays$RNA@data %>% as.data.frame() 
stria.data$gene <- rownames(stria.data)
stria.data.long <- stria.data %>%
  gather(cell,exp,-gene) %>% 
  left_join(stria.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(stria.data.long, file.path(Stria, "stria.data.long.tsv"), sep="\t", col.names=T)


##--------------- stria_im.rds  
stria_im <- readRDS(file.path(Stria, "stria_im.rds"))

#- cell ID to cell type
stria_im.meta <- stria_im@meta.data %>% select(celltype)
stria_im.meta$celltype <- gsub(" ","_", stria_im.meta$celltype)
stria_im.meta$cell <- rownames(stria_im.meta)

#- count
stria_im.counts <- stria_im@assays$RNA@counts %>% as.data.frame()
stria_im.counts$gene <- rownames(stria_im.counts)
stria_im.counts.long <- stria_im.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(stria_im.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(stria_im.counts.long, file.path(Stria, "stria_im.counts.long.tsv"), sep="\t", col.names=T)

#- normalized expression
stria_im.data <- stria_im@assays$RNA@data %>% as.data.frame() 
stria_im.data$gene <- rownames(stria_im.data)
stria_im.data.long <- stria_im.data %>%
  gather(cell,exp,-gene) %>% 
  left_join(stria_im.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(stria_im.data.long, file.path(Stria, "stria_im.data.long.tsv"), sep="\t", col.names=T)


##--------------- sgn.rds  
SGN <- "/home/wulab/scRNA/scRNA/lu"
# SGN <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_SGN"
sgn <- readRDS(file.path(SGN, "sgn.rds"))

#- cell ID to cell type
sgn.meta <- sgn@meta.data %>% select(celltype)
sgn.meta$celltype <- gsub(" ","_", sgn.meta$celltype)
sgn.meta$cell <- rownames(sgn.meta)

#- count
sgn.counts <- sgn@assays$RNA@counts %>% as.data.frame()
sgn.counts$gene <- rownames(sgn.counts)
sgn.counts.long <- sgn.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(sgn.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(sgn.counts.long, file.path(SGN, "sgn.counts.long.tsv"), sep="\t", col.names=T)

#- normalized expression
sgn.data <- sgn@assays$RNA@data %>% as.data.frame() 
sgn.data$gene <- rownames(sgn.data)
sgn.data.long <- sgn.data %>%
  gather(cell,exp,-gene) %>% 
  left_join(sgn.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(sgn.data.long, file.path(SGN, "sgn.data.long.tsv"), sep="\t", col.names=T)


##--------------- sgn_im.rds  
sgn_im <- readRDS(file.path(SGN, "sgn_im.rds"))

#- cell ID to cell type
sgn_im.meta <- sgn_im@meta.data %>% select(celltype)
sgn_im.meta$celltype <- gsub(" ","_", sgn_im.meta$celltype)
sgn_im.meta$cell <- rownames(sgn_im.meta)

#- count
sgn_im.counts <- sgn_im@assays$RNA@counts %>% as.data.frame()
sgn_im.counts$gene <- rownames(sgn_im.counts)
sgn_im.counts.long <- sgn_im.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(sgn_im.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(sgn_im.counts.long, file.path(SGN, "sgn_im.counts.long.tsv"), sep="\t", col.names=T)

#- normalized expression
sgn_im.data <- sgn_im@assays$RNA@data %>% as.data.frame() 
sgn_im.data$gene <- rownames(sgn_im.data)
sgn_im.data.long <- sgn_im.data %>%
  gather(cell,exp,-gene) %>% 
  left_join(sgn_im.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))
# output
fwrite(sgn_im.data.long, file.path(SGN, "sgn_im.data.long.tsv"), sep="\t", col.names=T)


# 2. Merge and Normalize count across cell types  
Stria <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_Stria"
SGN <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_SGN"

#- read in the exp.matrix: gene_x_cellTypes
sgn.counts.long      <- fread(file.path(SGN, "sgn.counts.long.tsv"), stringsAsFactors=F, data.table=F)
sgn_im.counts.long   <- fread(file.path(SGN, "sgn_im.counts.long.tsv"), stringsAsFactors=F, data.table=F)
stria.counts.long    <- fread(file.path(Stria, "stria.counts.long.tsv"), stringsAsFactors=F, data.table=F)
stria_im.counts.long <- fread(file.path(Stria, "stria_im.counts.long.tsv"), stringsAsFactors=F, data.table=F)

#- merge
dat <- rbind(sgn.counts.long,
             sgn_im.counts.long,
             stria.counts.long,
             stria_im.counts.long) %>% 
  select(-exp_tpm) %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) 

#- keep only mouse to human 1to1 mapped orthologs
dir <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision"
m2h <- fread(file.path(dir, "mouse_human_homologs.txt"))
names(m2h) <- c("MumSYM", "HumSYM")  

# dat <- dat %>% filter(gene %in% m2h$MumSYB) 
dat <- merge(dat, m2h, by.x="gene", by.y="MumSYM")

#- fill in NAs with 0:
tmp <- dat %>% pivot_wider(id_cols=c(gene, HumSYM), names_from=celltype, values_from=exp, values_fill=0, values_fn=list(exp=sum))
dat <- tmp %>% pivot_longer(cols=-c(gene, HumSYM), names_to="celltype", values_to="exp")

#- remove duplicated genes
tmp.dup <- dat %>% add_count(gene) 

#- remove genes not expressed in any cell type
tmp.noexp <- dat %>% 
  group_by(gene) %>% 
  summarise(sum_exp=sum(exp)) %>%
  filter(sum_exp==0)

dat <- dat %>% filter(!gene %in% tmp.noexp$gene)

#- add up count for duplicated cell types
dat <- dat %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp)) %>%
  ungroup()

# 3. Calculate specificity  
dat <- dat %>%
  group_by(gene) %>%
  mutate(specificity=exp_tpm/sum(exp_tpm)) %>%
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
# Filter out genes with expression below 1 TPM.
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon")
if (file.exists("MAGMA/top10.txt")) {file.remove("MAGMA/top10.txt")}
dat %>% filter(exp_tpm>1) %>% magma_top10("celltype")

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon")
dat %>% filter(exp_tpm>1) %>% ldsc_bedfile("celltype")
control <- as.data.frame(unique(dat$HumSYM))
readr::write_tsv(control,  "LDSC/control.bed", col_names=F)


