#- Trpchevska et al.
library(data.table)

setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/UKB/hearing-loss")
df <- fread("hl_ma_nov18_v2_summarystats_09122021.txt.gz")

# rename columns for mrmega in format
old_col <- c("BP", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P.value")
new_col <- c("POS", "A1", "A2", "freq", "beta", "SE", "p")
setnames(df, old_col, new_col)

is_sci <- grepl("[eE]", as.character(df$POS))
df[is_sci, ] 
df$POS <- as.integer(round(df$POS))

out_col <- df[, .(SNP, CHR, POS, A1, A2, freq, beta, SE, p, N)]
fwrite(out_col, file="AJHG_EUR_reformat.gz", sep="\t", row.names=FALSE)

