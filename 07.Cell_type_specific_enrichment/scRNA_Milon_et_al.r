#============================================================#
#
#  scRNA expression treatment of Milon
#
#============================================================#
#----------------------- stria treat -------------------------/
library(Seurat)
library(dplyr)
library(Matrix)
library(future)

set.seed(1234)
plan(sequential)
options(future.globals.maxSize = 16 * 1024^3)

Stria <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_Stria"
sample_dirs <- list.dirs(Stria, full.names=TRUE, recursive=FALSE)
objs_list <- list()
for (i in seq_along(sample_dirs)) {
  sample_path <- sample_dirs[i]
  sample_name <- basename(sample_path)
  sample_cond <- sub(".*_(naive|noise)\\d+$", "\\1", sample_name)
  message("Reading sample: ", sample_name)
  
  if (file.exists(file.path(sample_path, "genes.tsv.gz")) && !file.exists(file.path(sample_path, "features.tsv.gz"))) {
    file.copy(file.path(sample_path, "genes.tsv.gz"), file.path(sample_path, "features.tsv.gz"))
  }
  
  counts <- Read10X(data.dir = sample_path)
  seu <- CreateSeuratObject(counts=counts, project=sample_name, min.cells=3, min.features=100)
  seu$sample_id <- paste0("sample", i)
  seu$condition <- sample_cond
  
  seu <- RenameCells(seu, add.cell.id=sample_name)
  objs_list[[sample_name]] = seu
}

stria <- Reduce(function(x, y) merge(x, y), objs_list)
saveRDS(stria, file.path(Stria, "stria_raw_merged.rds"))

# ---------- 2) QC ----------
# stria <- readRDS(file.path(Stria, "stria_raw_merged.rds"))
DefaultAssay(stria) <- "RNA"
stria[["percent.mt"]] <- PercentageFeatureSet(stria, pattern = "(?i)^mt-")
stria[["percent.hb"]] <- PercentageFeatureSet(stria,pattern = "^Hba\\-a1$|^Hba\\-a2$|^Hbb\\-bh1$|^Hbb\\-bs$|^Hbb\\-bt$",assay   = "RNA")
stria <- subset(stria, subset=percent.hb < 1 & percent.mt < 25 & nFeature_RNA >= 1000 & nFeature_RNA <= 6000)
keep_features <- rowSums(GetAssayData(stria, slot = "counts") > 0) >= 20
stria <- subset(stria, features = rownames(stria)[keep_features])
saveRDS(stria, file.path(Stria, "stria_qc_filter.rds"))

# ---------- 3) sctransform → PCA/UMAP/cluster/annotation ----------
stria <- SCTransform(stria, verbose=FALSE, conserve.memory=TRUE, method = if (requireNamespace("glmGamPoi", quietly=TRUE)) "glmGamPoi" else "poisson")
stria <- RunPCA(stria, verbose=FALSE)
stria <- RunUMAP(stria, dims=1:25, verbose=FALSE)
stria <- FindNeighbors(stria, dims=1:25)
stria <- FindClusters(stria, resolution=0.6)

DimPlot(stria, label=TRUE)
Marginal_cells       = c("Abcg1","Heyl","Kcne1","Kcnq1")
Intermediate_Cells   = c("Cd44","Kcnj13","Met","Nrp2")
Basal_Cells          = c("Cldn11","Nr2f2","Sox8","Tjp1")
Fibrocytes           = c("Car3","Coch","Gm525","Igfbp2")
# Spindle/Root Cells
Spindle_Root_Cells   = c("Cldn9","Kcnj16","P2rx2","Slc26a4")

FeaturePlot(stria, Marginal_cells) # 1 6 27
DotPlot(object = stria, features = Marginal_cells)
FeaturePlot(stria, Intermediate_Cells)  # 0 2 4 5 9 12
DotPlot(object = stria, features = Intermediate_Cells)
FeaturePlot(stria, Basal_Cells) # 7 
DotPlot(object = stria, features = Basal_Cells)
FeaturePlot(stria, Fibrocytes) # 8 19 25
DotPlot(object = stria, features = Fibrocytes)
FeaturePlot(stria, Spindle_Root_Cells) # 14 15
DotPlot(object = stria, features = Spindle_Root_Cells)

# im
Monocytes            = c("Adgre1","Cd14","Cd68","Cx3cr1")
Neutrophils          = c("Lcn2","Ly6g")
B_Cells              = c("Cd19","Cd79a","Cd79b")

FeaturePlot(stria, Monocytes) # 20
DotPlot(object = stria, features = Monocytes)
FeaturePlot(stria, Neutrophils)  # 23
DotPlot(object = stria, features = Neutrophils)
FeaturePlot(stria, B_Cells) # 23
DotPlot(object = stria, features = B_Cells)

stria$celltype = plyr::mapvalues(stria$seurat_clusters,
                                  from = 0:27,
                                  to=c(
                                    "Intermediate_Cells", # 0
                                    "Marginal_cells",     # 1
                                    "Intermediate_Cells", # 2
                                    "unknow",             # 3
                                    "Intermediate_Cells", # 4
                                    "Intermediate_Cells", # 5
                                    "Marginal_cells",     # 6
                                    "Basal_Cells",        # 7
                                    "Fibrocytes",         # 8
                                    "Intermediate_Cells", # 9
                                    "unknow",             # 10
                                    "unknow",             # 11 -
                                    "Intermediate_Cells", # 12
                                    "unknow",             # 13 -
                                    "Spindle_Root_Cells", # 14
                                    "Spindle_Root_Cells", # 15
                                    "unknow",             # 16
                                    "unknow",             # 17
                                    "unknow",             # 18 -
                                    "Fibrocytes",         # 19
                                    "Monocytes",          # 20
                                    "unknow",             # 21
                                    "unknow",             # 22 -
                                    "imm",                # 23
                                    "unknow",             # 24
                                    "Fibrocytes",         # 25
                                    "unknow",             # 26
                                    "Marginal_cells"))    # 27


# ---------- 4) extract 8 target cell types，rerun PCA/UMAP/cluster ----------
five_cell <- c("Marginal_cells", "Intermediate_Cells", "Basal_Cells", "Fibrocytes","Spindle_Root_Cells")
imm_cell <- c("Monocytes", "imm")
stria5 <- stria[, stria@meta.data$celltype %in% five_cell]
stria5 <- RunPCA(stria5, verbose=FALSE)
stria5 <- RunUMAP(stria5, dims=1:25, verbose=FALSE)
stria5 <- FindNeighbors(stria5, dims=1:25)
stria5 <- FindClusters(stria5, resolution=0.6)
stria5@meta.data$celltype <- droplevels(stria5@meta.data$celltype)
saveRDS(stria5, file.path(Stria, "stria.rds"))

# immu
stria_im <- stria[, stria@meta.data$celltype %in% imm_cell]
stria_im <- RunPCA(stria_im, verbose=FALSE)
stria_im <- RunUMAP(stria_im, dims=1:25, verbose=FALSE)
stria_im <- FindNeighbors(stria_im, dims=1:25)
stria_im <- FindClusters(stria_im, resolution=0.6)
DimPlot(stria_im, label=TRUE)

Monocytes            = c("Adgre1","Cd14","Cd68","Cx3cr1")
Neutrophils          = c("Lcn2","Ly6g")
B_Cells              = c("Cd19","Cd79a","Cd79b")

FeaturePlot(stria_im, Monocytes) # 0 1 3 4
DotPlot(object = stria_im, features = Monocytes)
FeaturePlot(stria_im, Neutrophils)  # 6
DotPlot(object = stria_im, features = Neutrophils)
FeaturePlot(stria_im, B_Cells) # 2
DotPlot(object = stria_im, features = B_Cells)

stria_im$celltype = plyr::mapvalues(stria_im$seurat_clusters,
                                  from = 0:6,
                                  to=c(
                                    "Monocytes", # 0
                                    "Monocytes",     # 1
                                    "B_Cells", # 2
                                    "Monocytes",   # 3
                                    "Monocytes", # 4
                                    "unknow", # 5
                                    "Neutrophils"))  # 6
                                    

imm_cell <- c("Monocytes", "B_Cells", "Neutrophils")
stria_imm <- stria_im[, stria_im@meta.data$celltype %in% imm_cell]
stria_imm <- RunPCA(stria_imm, verbose=FALSE)
stria_imm <- RunUMAP(stria_imm, dims=1:25, verbose=FALSE)
stria_imm <- FindNeighbors(stria_imm, dims=1:25)
stria_imm <- FindClusters(stria_imm, resolution=0.6)
stria_imm@meta.data$celltype <- droplevels(stria_imm@meta.data$celltype)
saveRDS(stria_imm, file.path(Stria, "stria_im.rds"))


#----------------------- SGN treat -------------------------/
library(Seurat)
library(dplyr)
library(Matrix)
library(future)

set.seed(1234)
plan(sequential)
options(future.globals.maxSize = 16 * 1024^3)

SGN <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_SGN"
sample_dirs <- list.dirs(SGN, full.names=TRUE, recursive=FALSE)
objs_list <- list()
for (i in seq_along(sample_dirs)) {
  sample_path <- sample_dirs[i]
  sample_name <- basename(sample_path)
  sample_cond <- sub(".*_(naive|noise)\\d+$", "\\1", sample_name)
  message("Reading sample: ", sample_name)
  
  if (file.exists(file.path(sample_path, "genes.tsv.gz")) && !file.exists(file.path(sample_path, "features.tsv.gz"))) {
    file.copy(file.path(sample_path, "genes.tsv.gz"), file.path(sample_path, "features.tsv.gz"))
  }
  
  counts <- Read10X(data.dir = sample_path)
  seu <- CreateSeuratObject(counts=counts, project=sample_name, min.cells=3, min.features=100)
  seu$sample_id <- paste0("sample", i)
  seu$condition <- sample_cond
  
  seu <- RenameCells(seu, add.cell.id=sample_name)
  objs_list[[sample_name]] = seu
}

sgn <- Reduce(function(x, y) merge(x, y), objs_list)
saveRDS(sgn, file.path(SGN, "sgn_raw_merged.rds"))

# ---------- 2) QC ----------
# sgn <- readRDS(file.path(SGN, "sgn_raw_merged.rds"))
DefaultAssay(sgn) <- "RNA"
sgn[["percent.mt"]] <- PercentageFeatureSet(sgn, pattern = "(?i)^mt-")
sgn[["percent.hb"]] <- PercentageFeatureSet(sgn, pattern = "^Hba\\-a1$|^Hba\\-a2$|^Hbb\\-bh1$|^Hbb\\-bs$|^Hbb\\-bt$", assay = "RNA")
sgn <- subset(sgn, subset=percent.hb < 1 & percent.mt < 25 & nFeature_RNA >= 1000 & nFeature_RNA <= 6000)
keep_features <- rowSums(GetAssayData(sgn, slot = "counts") > 0) >= 20
sgn <- subset(sgn, features = rownames(sgn)[keep_features])
saveRDS(sgn, file.path(SGN, "sgn_qc_filter.rds"))

# ---------- 3) sctransform → PCA/UMAP/cluster/annotation ----------
sgn <- SCTransform(sgn, verbose=FALSE, conserve.memory=TRUE, method = if(requireNamespace("glmGamPoi", quietly=TRUE)) "glmGamPoi" else "poisson")
sgn <- RunPCA(sgn, verbose=FALSE)
sgn <- RunUMAP(sgn, dims=1:25, verbose=FALSE)
sgn <- FindNeighbors(sgn, dims=1:25)
sgn <- FindClusters(sgn, resolution=0.6)

DimPlot(sgn, label=TRUE)

Schwann_Cells        = c("Mbp","Mpz","Mpzl1","Pmp22")
Type_1               = c("Chgb","Kcnc3","Nefl","Scn4b", "Tubb3")
Type_2               = c("Gata3","Mafb","Ngfr","Prph", "Th")
Monocytes            = c("Cd68","Cx3cr1")
Neutrophils          = c("Lcn2","Ly6g")

FeaturePlot(sgn, Schwann_Cells) # 3 29
DotPlot(object = sgn, features = Schwann_Cells)
FeaturePlot(sgn, Type_1)  # 0 6 15 
DotPlot(object = sgn, features = Type_1)
FeaturePlot(sgn, Type_2) # 12
DotPlot(object = sgn, features = Type_2)

# im
FeaturePlot(sgn, Monocytes) # 13
DotPlot(object = sgn, features = Monocytes)
FeaturePlot(sgn, Neutrophils) # 19
DotPlot(object = sgn, features = Neutrophils)

sgn$celltype = plyr::mapvalues(sgn$seurat_clusters,
                                    from = 0:29,
                                    to=c(
                                      "Type_1",        # 0
                                      "unknow",        # 1
                                      "unknow",        # 2
                                      "Schwann_Cells", # 3
                                      "unknow",        # 4
                                      "unknow",        # 5
                                      "Type_1",        # 6
                                      "unknow",        # 7
                                      "unknow",        # 8
                                      "unknow",        # 9
                                      "unknow",        # 10
                                      "unknow",        # 11
                                      "Type_2",        # 12
                                      "Monocytes",     # 13
                                      "unknow",        # 14
                                      "Type_1",        # 15
                                      "unknow",        # 16
                                      "unknow",        # 17
                                      "unknow",        # 18
                                      "Neutrophils",   # 19
                                      "unknow",        # 20
                                      "unknow",        # 21
                                      "unknow",        # 22
                                      "unknow",        # 23
                                      "unknow",        # 24
                                      "unknow",        # 25
                                      "unknow",        # 26
                                      "Type_1",        # 27
                                      "unknow",        # 28
                                      "Schwann_Cells"  # 29
                                      ))  

# ---------- 4) extract SGN and Schwann cells，rerun PCA/UMAP/cluster/annotation  ----------
SGN_Schwann <- c("Type_1", "Type_2", "Schwann_Cells")
sgn_schwann <- sgn[, sgn@meta.data$celltype %in% SGN_Schwann]
sgn_schwann <- RunPCA(sgn_schwann, verbose=FALSE)
sgn_schwann <- RunUMAP(sgn_schwann, dims=1:25, verbose=FALSE)
sgn_schwann <- FindNeighbors(sgn_schwann, dims=1:25)
sgn_schwann <- FindClusters(sgn_schwann, resolution=0.6)
sgn_schwann@meta.data$celltype <- droplevels(sgn_schwann@meta.data$celltype)
DimPlot(sgn_schwann, label=TRUE)

FeaturePlot(sgn_schwann, Schwann_Cells) # 0 4 12
DotPlot(object = sgn_schwann, features = Schwann_Cells)
FeaturePlot(sgn_schwann, Type_2) # 6 10
DotPlot(object = sgn_schwann, features = Type_2)

Type_1A              = c("B3gat1","Calb2","Obscn","Pcdh20")
Type_1B              = c("Calb1","Runx1","Ttn")
Type_1C              = c("Grm8","Hmcn1","Kcnip2", "Lypd1", "Pou4f1")
FeaturePlot(sgn_schwann, Type_1A) # 1 2 5 7
DotPlot(object = sgn_schwann, features = Type_1A)
FeaturePlot(sgn_schwann, Type_1B)  # 3 8 11
DotPlot(object = sgn_schwann, features = Type_1B)
FeaturePlot(sgn_schwann, Type_1C) # 9
DotPlot(object = sgn_schwann, features = Type_1C)

sgn_schwann$celltype = plyr::mapvalues(sgn_schwann$seurat_clusters,
                                from = 0:12,
                                to=c(
                                  "Schwann_Cells", # 0
                                  "Type_1A",       # 1
                                  "Type_1A",       # 2
                                  "Type_1B",       # 3
                                  "Schwann_Cells", # 4
                                  "Type_1A",       # 5
                                  "Type_2",        # 6
                                  "Type_1A",       # 7
                                  "Type_1B",       # 8
                                  "Type_1C",       # 9
                                  "Type_2",        # 10
                                  "Type_1B",       # 11
                                  "Schwann_Cells"  # 12
                                ))  
saveRDS(sgn_schwann, file.path(SGN, "sgn.rds"))


# im
SGN_im<- c("Monocytes", "Neutrophils")
sgn_im <- sgn[, sgn@meta.data$celltype %in% SGN_im]
sgn_im <- RunPCA(sgn_im, verbose=FALSE)
sgn_im <- RunUMAP(sgn_im, dims=1:25, verbose=FALSE)
sgn_im <- FindNeighbors(sgn_im, dims=1:25)
sgn_im <- FindClusters(sgn_im, resolution=0.6)
sgn_im@meta.data$celltype <- droplevels(sgn_im@meta.data$celltype)

DimPlot(sgn_im, label=TRUE)
saveRDS(sgn_im, file.path(SGN, "sgn_im.rds"))

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


#### scDRS file 
library(Seurat)
library(stringr)
library(SeuratDisk) 

Stria <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_Stria"
SGN <- "F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon_et_al_GSE168041_RAW_SGN"

stria <- file.path(Stria, "stria.rds")
stria_im <- file.path(Stria, "stria_im.rds")
sgn <- file.path(SGN, "sgn.rds")
sgn_im <- file.path(SGN, "sgn_im.rds")

rds_files <- c(stria, stria_im, sgn, sgn_im)

objs <- lapply(seq_along(rds_files), function(i){
  obj <- readRDS(rds_files[i])
  DefaultAssay(obj) <- "RNA"

  obj
  })

combined <- merge(x = objs[[1]], y = objs[-1])

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, normalization.method="LogNormalize", verbose=FALSE)

# delete dup genes
if (any(duplicated(rownames(combined)))) {
  rownames(combined) <- make.unique(rownames(combined))
}

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/scDT/Milon")
SaveH5Seurat(combined, filename = "Milon_et_al.h5seurat", overwrite = TRUE)
Convert("Milon_et_al.h5seurat", dest = "h5ad", assay = "RNA", overwrite = TRUE)


# check in scanpy
import scanpy as sc, numpy as np

ad = sc.read_h5ad("Milon_et_al.h5ad")

if "counts" not in ad.layers.keys():
    if ad.raw is not None:
        ad.layers["counts"] = ad.raw.X.copy()
    else:
        ad.layers["counts"] = ad.X.copy()

if "n_counts" not in ad.obs:
    ad.obs["n_counts"] = np.asarray(ad.layers["counts"].sum(axis=1)).ravel()
if "n_genes" not in ad.obs:
    ad.obs["n_genes"] = np.asarray((ad.layers["counts"]>0).sum(axis=1)).ravel()

ad.write_h5ad("Milon_et_al_ready.h5ad", compression="gzip")