library(Seurat)
library(tidyverse)
library(dplyr)

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Xu_et_al_P28_GSE202920_RAW")

#P28_1#####
counts2mp28_1 <- Read10X("p28_1/")
scrna_2mp28_1 <- CreateSeuratObject(counts2mp28_1, min.cells = 3)
scrna_2mp28_1@meta.data[["orig.ident"]] = "p28_1"

scrna_2mp28_1[["percent.mt"]] <- PercentageFeatureSet(scrna_2mp28_1, pattern = "^mt-")
scrna_2mp28_1[["percent_Rp.sl."]] <- PercentageFeatureSet(scrna_2mp28_1, pattern = "^Rp[sl]")

HB.genes <- c("Hba-a1","Hba-a2","Hbb-bh1","Hbb-bs","Hbb-bt")
HB_m <- match(HB.genes, rownames(scrna_2mp28_1@assays$RNA)) 
Hb.genes <- rownames(scrna_2mp28_1@assays$RNA)[HB_m] 
Hb.genes <- Hb.genes[!is.na(Hb.genes)] 
scrna_2mp28_1[["percent_Hblist"]]<-PercentageFeatureSet(scrna_2mp28_1, features=Hb.genes)

ScRNA_2mp28_1 <- subset(scrna_2mp28_1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &nCount_RNA < 15000 & percent.mt < 15)


#p28_2#####
counts2mp28_2 <- Read10X("p28_2/")
scrna_2mp28_2 <- CreateSeuratObject(counts2mp28_2, min.cells = 3)
scrna_2mp28_2@meta.data[["orig.ident"]] = "p28_2"

scrna_2mp28_2[["percent.mt"]] <- PercentageFeatureSet(scrna_2mp28_2, pattern = "^mt-")
scrna_2mp28_2[["percent_Rp.sl."]] <- PercentageFeatureSet(scrna_2mp28_2, pattern = "^Rp[sl]")

HB.genes <- c("Hba-a1","Hba-a2","Hbb-bh1","Hbb-bs","Hbb-bt")
HB_m <- match(HB.genes, rownames(scrna_2mp28_2@assays$RNA)) 
Hb.genes <- rownames(scrna_2mp28_2@assays$RNA)[HB_m] 
Hb.genes <- Hb.genes[!is.na(Hb.genes)] 
scrna_2mp28_2[["percent_Hblist"]]<-PercentageFeatureSet(scrna_2mp28_2, features=Hb.genes)

ScRNA_2mp28_2 <- subset(scrna_2mp28_2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 15000 & percent.mt < 15)


# P28 Seurat
ScRNA_p28 <- list()
ScRNA_p28 <-  list (ScRNA_2mp28_1,ScRNA_2mp28_2)

rm(scrna_2mp28, scrna_2mp28_1,scrna_2mp28_2,
   ScRNA_2mp28,ScRNA_2mp28_1,ScRNA_2mp28_2,
   counts2mp28,counts2mp28_1,counts2mp28_2)
rm(feats, HB_m,HB.genes,Hb.genes)
gc()

#Integration####
ScRNA_p28 <- lapply(X = ScRNA_p28, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ScRNA_p28)

ScRNA_p28.anchors <- FindIntegrationAnchors(object.list = ScRNA_p28, anchor.features = features)
ScRNA_p28 <- IntegrateData(anchorset = ScRNA_p28.anchors)
table(ScRNA_p28@meta.data$orig.ident)

save(ScRNA_p28, file = "Integrated_ScRNA_p28")
#Integrated analysis####
DefaultAssay(ScRNA_p28) <- "integrated"
ScRNA_p28 = ScaleData(ScRNA_p28)
ScRNA_p28 <- RunPCA(ScRNA_p28, npcs = 30, verbose = T)

DimPlot(ScRNA_p28, reduction = "pca", group.by="orig.ident") 
ElbowPlot(ScRNA_p28, ndims=30, reduction="pca") 

pc.num=1:30
ScRNA_p28 <- FindNeighbors(ScRNA_p28, reduction = "pca", dims = pc.num)
ScRNA_p28 <- FindClusters(ScRNA_p28, resolution = 0.8)
ScRNA_p28 <- RunUMAP(ScRNA_p28, reduction = "pca", dims = pc.num)
ScRNA_p28 <- RunTSNE(ScRNA_p28, dims = pc.num)

ScRNA_p28$celltype = plyr::mapvalues(ScRNA_p28$seurat_clusters,
                                     from = 0:27,
                                     to=c(
                                       "BaC",           # 0
                                       "TBC",           # 1
                                       "ISC",           # 2
                                       "TBC",           # 3
                                       "RM",            # 4
                                       "MP",            # 5
                                       "HeC",           # 6
                                       "IPhC_IB",       # 7
                                       "FC_Abcc9+",     # 8
                                       "BaC",           # 9
                                       "OHC",           # 10
                                       "OSC",           # 11
                                       "SPC_RC",        # 12
                                       "FC_Coch+_Spp1+",# 13
                                       "RBC",           # 14
                                       "FC_Anxa1",      # 15
                                       "BaC",           # 16
                                       "SGN",           # 17
                                       "HeC",           # 18
                                       "IC",            # 19
                                       "SC_SGC",        # 20
                                       "DC_PC",         # 21
                                       "SPC_RC",        # 22
                                       "IHC",           # 23
                                       "BC",            # 24
                                       "SGN",           # 25
                                       "FC_Bglap+",     # 26
                                       "MC"             # 27
                                     ))  

saveRDS(ScRNA_p28, "ScRNA_p28.rds")
### save h5ad for scDRS 
library(SeuratDisk) 

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Xu")
obj <- readRDS("ScRNA_p28.rds")
ScRNA_p28$celltype <- as.character(ScRNA_p28$celltype)   

SaveH5Seurat(obj, filename = "Xu_et_al.h5seurat", overwrite = TRUE)
Convert("Xu_et_al.h5seurat", dest = "h5ad", assay = "RNA", overwrite = TRUE)


### for LDSC-SEG MAGMA
p28.meta <- ScRNA_p28@meta.data %>% select(celltype)
p28.meta$cell <- rownames(p28.meta)

#- count
p28.counts <- ScRNA_p28@assays$RNA@counts %>% as.data.frame()
p28.counts$gene <- rownames(p28.counts)
p28.counts.long <- p28.counts %>%
  gather(cell,exp,-gene) %>% 
  left_join(p28.meta, by="cell") %>%
  group_by(gene, celltype) %>%
  summarise(exp=sum(exp)) %>%
  group_by(celltype) %>%
  mutate(exp_tpm=exp*1e6/sum(exp))


dat <- p28.counts.long
#- keep only mouse to human 1to1 mapped orthologs
dir <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision"
m2h <- fread(file.path(dir, "mouse_human_homologs.txt"))
names(m2h) <- c("MumSYM", "HumSYM")  

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
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Xu")
if (file.exists("MAGMA/top10.txt")) {file.remove("MAGMA/top10.txt")}
dat %>% filter(exp_tpm>1) %>% magma_top10("celltype")

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Xu")
dat %>% filter(exp_tpm>1) %>% ldsc_bedfile("celltype")
control <- as.data.frame(unique(dat$HumSYM))
readr::write_tsv(control,  "LDSC/control.bed", col_names=F)
