library(data.table)
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss/hum0197.v3.BBJ.HL.v1")
df <- fread("GWASsummary_Hearing_Loss_Japanese_SakaueKanai2020.auto.txt.gz")

# rename columns for mrmega in format
old_col <- c("SNPID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value")
new_col <- c("SNP", "A2", "A1", "freq", "beta", "SE", "p")
setnames(df, old_col, new_col)

is_sci <- grepl("[eE]", as.character(df$POS))
df[is_sci, ] 
df$POS <- as.integer(round(df$POS))

out_col <- df[, .(SNP, CHR, POS, A1, A2, freq, beta, SE, p, N)]

out_col = out_col[out_col$CHR != "X", ]
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss")
fwrite(out_col, file="BBJ_EAS_reformat.gz", sep="\t", row.names=FALSE)

