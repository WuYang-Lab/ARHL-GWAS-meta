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


