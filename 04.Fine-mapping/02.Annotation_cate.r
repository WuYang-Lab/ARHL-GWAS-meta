#=============================================================#
# # # generate SNP annotation of  Number of 1-SNP credible sets
#=============================================================#
# read local credible set
setwd("/public/share/wchirdzhq2022/Wulab_share/sll/ARHL/NC_sup_test/04.gctb/gwfm")
library(data.table)
library(dplyr)

lcs = fread('All_MVP_Trpchevska_De-Angelis_BBJ_SbayesRC_finemapping.lcs')[, .(Size, SNP)][, SNP := gsub("[\t\r\n]", "", SNP)]; 
lcs_snp = lcs[Size==1, SNP]
baseline2 = fread('/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt')
lcs_snp_annot = baseline2[SNP %in% lcs_snp]
fwrite(lcs_snp_annot, 'finemapping_snp.txt', col.names=TRUE, sep='\t')
