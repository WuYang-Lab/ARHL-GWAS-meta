############################################
##   cross ancestry meta of hearing loss  ##
############################################

args       <- commandArgs(TRUE)
infile     <- args[1]
outfile    <- args[2]

library(dplyr)
library(stringr)
library(data.table)
df <- fread(infile)

# rename columns for mrmega in format
old_col <- c("SNP_ID", "chrom", "pos", "ref", "alt", "af", "num_samples", "or", "pval")
new_col <- c("SNP", "CHR", "POS", "A2", "A1", "freq", "N", "OR", "p")
setnames(df, old_col, new_col)

# OR to BETA and CI to SE
df$beta <- log(df$OR)

# if beta=0, then se=mean(allse)
se_mean <- df %>%
  filter(beta!=0) %>%
  mutate(se = sqrt(beta^2 / qchisq(p, 1, lower.tail=FALSE))) %>%
  summarise(mean_se = mean(se, na.rm = TRUE)) %>%
  pull(mean_se)
df$SE <- ifelse(df$beta==0, se_mean, sqrt(df$beta^2/qchisq(df$p, 1, lower.tail=F)))
df$Idx <- 1
# save the file
out_col <- df[, .(SNP, CHR, POS, A1, A2, freq, beta, SE, p, N, Idx)]
fwrite(out_col, file=paste0(outfile, ".gz"), sep="\t", row.names=FALSE)

