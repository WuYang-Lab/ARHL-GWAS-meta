setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta")
infile="ARHL_MVP_AJHG_BBJ_reformatMETAL.gz"
plot_mhn(infile, CHR="CHR", BP="POS", p="p")

# # # plot qq
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript
library(CMplot)
library(data.table)
plot_qq <- function(infile, SNP="SNP", CHR="CHR", BP="BP", P="P")

data <- fread("ARHL_MVP_AJHG_BBJ_reformatMETAL_addchr.gz")
pltdt <- subset(data, select=c(SNP, CHR, POS, p))
colnames(pltdt) <- c("SNP", "CHR", "BP", "P")
plt_filter <- pltdt[pltdt$P > 0, ]
plt_filter$P <- -log10(plt_filter$P)

CMplot(plt_filter, plot.type="q", LOG10=FALSE, file="jpg", file.name="ARHL_qqplot", dpi=500, 
      verbose=TRUE, cex=0.8, signal.cex=0.8, width=10, height=10, 
      threshold.col="red", threshold.lty=2, conf.int=TRUE, conf.int.col=NULL, file.output=TRUE)