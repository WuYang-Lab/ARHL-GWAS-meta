library(dplyr)
library(stringr)
library(data.table)
setwd("/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/01.meta")
df <- fread( "All_MVP_Trpchevska_De-Angelis_BBJ1.tbl")[, c("Chromosome", "Position", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize")]
names(df) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "HetISq", "HetChiSq", "HetDf", "HetPVal", "N")

df_rfm <- df %>%
    filter(HetDf > 0) %>%
    mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2)) %>%
    group_by(SNP) %>%
    slice_min(order_by = p, with_ties = FALSE) %>%
    ungroup()


out_col1 <- df_rfm[, c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")]
fwrite(out_col1, "All_MVP_Trpchevska_De-Angelis_BBJ_filter_chr.gz", sep="\t")

out_col2 <- df_rfm[, c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")]
fwrite(out_col2, "All_MVP_Trpchevska_De-Angelis_BBJ_filter.gz", sep="\t")
fwrite(out_col2, "All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt", sep="\t")
