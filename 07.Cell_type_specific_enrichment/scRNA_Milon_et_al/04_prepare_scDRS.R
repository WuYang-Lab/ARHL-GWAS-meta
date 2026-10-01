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


