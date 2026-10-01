library(dplyr)
library(stringr)
library(data.table)

setwd("/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta")
files <- c("EAS_MVP_BBJ", "EUR_MVP_Trpchevska_De-Angelis")

for (file in files){
  df <- fread(paste0(file, "1.tbl"))[, c("Chromosome", "Position", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "HetDf", "TotalSampleSize")]
  names(df) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "HetDf", "N")

  # reformat and filter the METAL out file
  df_rfm <- df %>%
    filter(HetDf > 0) %>%
    mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))

  out_col <- df_rfm[, c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")]
  fwrite(out_col, file=paste0(file, "_METALfilter.gz"), sep="\t", row.names=FALSE)
}


dt = fread("EUR_MVP_Trpchevska_De-Angelis_METALfilter.gz")
intercept = 1.1383 
sca = sqrt(intercept)
dt[, Z := beta / SE]
dt[, `:=`(SE_int = SE * sca, Z_int = beta / (SE * sca), P_int = 2 * pnorm(-abs(beta / (SE * sca))))]
out = dt[, .(CHR, POS, SNP, A1, A2, freq, beta, SE_int, P_int, N)]
names(out) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")
fwrite(out, "EUR_MVP_Trpchevska_De-Angelis_METALfilter_interceptAdj.gz", sep="\t")


