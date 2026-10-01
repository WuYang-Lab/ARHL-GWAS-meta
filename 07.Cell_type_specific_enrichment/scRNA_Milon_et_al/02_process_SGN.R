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

