############################################
##   cross ancestry meta of hearing loss  ##
############################################
#=======================================================#
#####  MVP hearing loss GWAS format convert         #####
#=======================================================#
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
#=======================================================
#! /bin/bash
RDIR="/public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript"
outdir="/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss"

for file in *.GIA.dbGaP.txt.gz; do
    id=$(echo "${file}" | cut -d '.' -f 4)
    ${RDIR} reformat_MVP.r ${file} ${outdir}/MVP_${id}_reformat
    qsubshcom "$cmd" 1 100G reformat 90:00:00 ""
done

#=======================================================#
###### BBJ hearing loss GWAS format convert        #####
#=======================================================#
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript

library(stringr)
library(data.table)

setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss/hum0197.v3.BBJ.HL.v1")
df <- fread("GWASsummary_Hearing_Loss_Japanese_SakaueKanai2020.auto.txt.gz")

# rename columns for mrmega in format
old_col <- c("SNPID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value")
new_col <- c("SNP", "A2", "A1", "freq", "beta", "SE", "p")
setnames(df, old_col, new_col)
df$Idx <- 1

out_col <- df[, .(SNP, CHR, POS, A1, A2, freq, beta, SE, p, N, Idx)]
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss")
fwrite(out_col, file="BBJ_EAS_reformat.gz", sep="\t", row.names=FALSE)
#=======================================================

#=======================================================#
#####  UKB hearing loss GWAS format convert         #####
#=======================================================#
#! /public/home/shilulu/Lu_mamba/Configs/envs/R4.2.0/bin/R

library(stringr)
library(data.table)

setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/UKB/hearing-loss")
df <- fread("hl_ma_nov18_v2_summarystats_09122021.txt.gz")

# rename columns for mrmega in format
old_col <- c("BP", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P.value")
new_col <- c("POS", "A1", "A2", "freq", "beta", "SE", "p")
setnames(df, old_col, new_col)
df$Idx <- 1

out_col <- df[, .(SNP, CHR, POS, A1, A2, freq, beta, SE, p, N, Idx)]
fwrite(out_col, file="AJHG_EUR_reformat.gz", sep="\t", row.names=FALSE)
#======================================================

#==========================================================#
#       *^  METAL for meta analysis  ^* 
#==========================================================#

#============ The EAS meta BBJ and MVP ============#
# METAL Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
TRACKPOSITIONS ON 
VERBOSE OFF
# Input columns:
CHROMOSOME CHR
MARKER SNP
ALLELE A1 A2
EFFECT beta
STDERR SE
FREQ freq
PVALUE p
WEIGHT N
# additional options
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
CUSTOMVARIABLE Indix 
LABEL Indix as Idx

PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss/BBJ_EAS_reformat.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_EAS_reformat.gz

OUTFILE ARHL_EAS_MVP_BBJ.test .tbl
ANALYZE HETEROGENEITY 
QUIT
#=========================
# run metal
METAL="/public/home/shilulu/software/metal/metal"
qsubshcom "$METAL EAS_MVP_BBJ_metal.conf" 1 100G METAL 90:00:00 ""
#=========================

#============ The EUR meta AJHG and MVP ============#
# METAL Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
# Input columns:
MARKER SNP
ALLELE A1 A2
EFFECT beta
STDERR SE
FREQ freq
PVALUE p
WEIGHT N
# additional options
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
CUSTOMVARIABLE Indix 
LABEL Indix as Idx

PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/UKB/hearing-loss/AJHG_EUR_reformat.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_EUR_reformat.gz

OUTFILE /public/home/shilulu/project_hearing-loss/new_run/ARHL_EUR_MVP_AJHG .tbl
ANALYZE HETEROGENEITY 
QUIT
#=========================
# run metal
METAL="/public/home/shilulu/software/metal/metal"
qsubshcom "$METAL EUR_MVP_AJHG_metal.conf" 1 100G METAL 90:00:00 ""
#=========================

#======================
# reformat and filter METAL result 
# the EAS and EUR, which filter the SNP only occured in one study
library(dplyr)
library(stringr)
library(data.table)

setwd("/public/home/shilulu/project_hearing-loss/new_run")
files <- c("ARHL_EAS_MVP_BBJ", "ARHL_EUR_MVP_AJHG")

METAL_rfm_flt <- function(file) {
  df <- fread(paste0(file, "1.tbl.gz"))
  out <- paste0(file, "_reformatMETAL")

  old_col = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize", "Indix")
  new_col = c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N", "Idx")
  setnames(df, old_col, new_col)
  # reformat the METAL out file
  df_rfm <- df %>%
    filter(Idx > 1) %>%
    mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2)) %>%
    mutate(Idx = 1)
  out_col <- df_rfm[, .(SNP, A1, A2, freq, beta, SE, p, N, Idx)]
  fwrite(out_col, file=paste0(out, ".gz"), sep="\t", row.names=FALSE)
}

for (file in files){
  METAL_rfm_flt(file)
}

#=============================#
#  meta for EAS EUR AFR AMR   #
# the mean result to analysis #
#=============================#
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
# Input columns:
MARKER SNP
ALLELE A1 A2
EFFECT beta
STDERR SE
FREQ freq
PVALUE p
WEIGHT N
# additional options
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
CUSTOMVARIABLE Indix 
LABEL Indix as Idx

PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_AFR_reformat.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_AMR_reformat.gz
PROCESS /public/home/shilulu/project_hearing-loss/new_run/ARHL_EAS_MVP_BBJ_reformatMETAL.gz
PROCESS /public/home/shilulu/project_hearing-loss/new_run/ARHL_EUR_MVP_AJHG_reformatMETAL.gz

OUTFILE /public/home/shilulu/project_hearing-loss/new_run/ARHL_MVP_AJHG_BBJ .tbl
ANALYZE HETEROGENEITY 
QUIT

#=========================
# run metal
METAL="/public/home/shilulu/software/metal/metal"
qsubshcom "$METAL MVP_AJHG_BBJ_metal.conf" 1 100G METAL 90:00:00 ""
gzip ARHL_MVP_AJHG_BBJ1.tbl
#==================================
# reformat and filter METAL result 
# and add chr pos for the result
#  rfm_flt_addchrpos_metal.r
#==================================
library(dplyr)
library(data.table)

# reformat the METAL out file
# filter the SNP which only occured at one study
METAL_rfm_flt <- function(file){
  df <- fread(paste0(file, "1.tbl.gz"))

  old_col = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize")
  new_col = c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")
  setnames(df, old_col, new_col)

  df_rfm <- df %>%
    filter(Indix > 1) %>%
    mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))
  out_col <- df_rfm[, .(SNP, A1, A2, freq, beta, SE, p, N)]
  return(out_col)
}

# add chr and pos information
add_chr_pos <- function(df) {
  setwd("/public/share/wchirdzhq2022/Wulab_share/dbSNP/GRCh37")
  db = fread(paste0(chr, ".txt"), header=FALSE, col.names=c("CHR", "POS", "SNP"))
  df_merge = merge(df, db, by="SNP", all=FALSE)
  return(df_merge)
}

setwd("/public/home/shilulu/project_hearing-loss/new_run")
metal_file = c("ARHL_MVP_AJHG_BBJ")

flt_metal = METAL_rfm_flt(metal_file)
fwrite(flt_metal, file=paste0(metal_file, "_reformatMETAL"), sep="\t", row.names=FALSE, quote=FALSE)

out_added = data.table()
for (chr in c(1:22)){
  out_added = rbind(out_added, add_chr_pos(flt_metal))
}
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta")
fwrite(out_added, file=paste0(metal_file, "_reformatMETAL_addchr.gz"), sep="\t", row.names=FALSE, quote=FALSE)
##########
qsubshcom "Rscript rfm_flt_addchrpos_metal.r" 1 100G METAL 90:00:00 ""
#############



###########################
## PLINK1.9  LD clumping ##
# hearing loss
###########################
# using eur reference genotype
eurdir="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur"
indir="/public/home/shilulu/project_hearing-loss/new_run/all_meta"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ld_clump"
cmd="plink \
    --clump ${indir}/ARHL_MVP_AJHG_BBJ_reformatMETAL.gz \
    --bfile ${eurdir}/g1000_eur \
    --clump-p1 5e-8 \
    --clump-p2 5e-8 \
    --clump-r2 0.05 \
    --clump-kb 1000 \
    --clump-snp-field SNP \
    --clump-field p \
    --out ${outdir}/ARHL_MVP_AJHG_BBJ_reformatMETAL"
qsubshcom "$cmd" 30 100G plink_clump 90:00:00 ""
