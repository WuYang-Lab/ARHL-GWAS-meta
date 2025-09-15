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

# hg38to19
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss")
library(data.table)

build_db <- function(dir = "/public/share/wchirdzhq2022/Wulab_share/dbSNP/GRCh37/") {
  db_all = rbindlist(lapply(1:22, function(chr) {
    f = file.path(dir, paste0(chr, ".txt"))
    fread(f, header = FALSE, col.names = c("CHR", "POS", "SNP"))
  }))
  return(db_all)
}
db_ref = build_db()

file38 = c("MVP_EAS_reformat.gz", "MVP_EUR_reformat.gz", "MVP_AFR_reformat.gz", "MVP_AMR_reformat.gz")

for (f in file38) {
  dt = fread(f)[, .(SNP, A1, A2, freq, beta, SE, p, N)]
  out_added = merge(dt, db_ref, by = "SNP", all = FALSE)
  setcolorder(out_added, c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N"))  
  out_name = sub("\\.gz$", "_hg19.gz", f)
  fwrite(out_added, out_name, sep = "\t")
}

#=======================================================#
###### BBJ hearing loss GWAS format convert        #####
#=======================================================#
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
#=======================================================


#=======================================================#
#####  EUR hearing loss GWAS format convert         #####
#=======================================================#
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

#- de-meta for De-Angelis et al.
gwasQC <- function(gwas, MAF=0.01){
    if (!requireNamespace("logger", quietly=TRUE)){ install.packages("logger") }
    lapply(c("dplyr","data.table","logger"), function(pkg) suppressWarnings(library(pkg,character.only=TRUE)))
    
    names(gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
    # beta != 0
    check_beta <- function(gwas){
        gwas$b = as.numeric(gwas$b)
        raw_n = nrow(gwas)
        gwas = gwas[!is.na(gwas$b) & is.finite(gwas$b) & gwas$b != 0, ]
        del_n = raw_n - nrow(gwas)

        return(list(n = del_n, gwas = gwas))
    }
    # 0 <= p <= 1
    check_p <- function(gwas){
        gwas$p = as.numeric(gwas$p)
        raw_n = nrow(gwas)
        gwas = gwas[!is.na(gwas$p) & is.finite(gwas$p) & gwas$p <= 1 & gwas$p >= 0, ]
        del_n = raw_n - nrow(gwas)

        return(list(n = del_n, gwas = gwas))
    }
    # MAF
    check_frq <- function(gwas, MAF){
        gwas$freq = as.numeric(gwas$freq)
        raw_n = nrow(gwas)
        gwas = gwas[!is.na(gwas$freq) & is.finite(gwas$freq) & gwas$freq > MAF & gwas$freq < 1-MAF, ]
        del_n = raw_n - nrow(gwas)

        return(list(n = del_n, gwas = gwas))
    }

    raw_n = nrow(gwas)
    log_info(sprintf("GWAS has %d SNPs before QC!", raw_n))
    log_info("Start to check SNP effect size")
    res = check_beta(gwas)
    gwas = res$gwas
    del_n = res$n 
    if(del_n > 0){
        log_info(sprintf("There were %d SNPs have deleted by abnormal effect size! (e.g., b=0 or b=NA)", del_n))
    } else {log_info("All SNPs effect size is normal!")}

    log_info("Start to check P value")
    res = check_p(gwas)
    gwas = res$gwas
    del_n = res$n
    if(del_n > 0){
        log_info(sprintf("There were %d SNPs have deleted by abnormal p value! (i.e., p < 0 or p > 1)", del_n))
    }else{ log_info("All SNPs p value is normal!") }

    log_info("Start to do MAF filter")
    res = check_frq(gwas, MAF)
    gwas = res$gwas
    del_n = res$n
    if(del_n > 0){
        log_info(sprintf("There were %d SNPs have deleted by MAF < %.4f! (i.e., %.4f < freq < %.4f)", del_n, MAF, MAF, 1 - MAF))
    } else {log_info(sprintf("All SNPs MAF > %4.f (i.e., %4.f < freq < %4.f)", MAF, MAF, 1 - MAF))}
    
    return(gwas)
}


de_meta <- function(gwas1, gwas2, diffreq=0.2){
  if (!requireNamespace("logger", quietly=TRUE)){ install.packages("logger") }
  lapply(c("dplyr","data.table","logger"), function(pkg) suppressWarnings(library(pkg,character.only=TRUE)))
  
  log_info("Meta summary is gwas1 and cohort want to delete is gwas2")
  log_info("Start to merge the gwas1 and gwas2 by SNP")
  gwas12 = merge(gwas1, gwas2, by="SNP", suffixes=c("_1", "_2"))

  raw_n = nrow(gwas12)
  log_info(sprintf("There were %d SNPs common in gwas1 and gwas2!", raw_n))
  gwas1Uniq = gwas1[gwas1$SNP %in% setdiff(gwas1$SNP, gwas2$SNP), ]
  log_info(sprintf("There were %d SNPs only occured in gwas1, will delete these SNPs!", nrow(gwas1Uniq)))

  log_info("Start to fix the effect allele in gwas1 and gwas2")
  gwas12[, b_2 := fifelse((A1_1==A1_2 & A2_1==A2_2), b_2, fifelse(A1_1==A2_2 & A2_1==A1_2, -b_2, NA_real_))]
  gwas12[, freq_2 := fifelse((A1_1==A1_2 & A2_1==A2_2), freq_2, fifelse(A1_1==A2_2 & A2_1==A1_2, 1-freq_2,NA_real_))]
  
  log_info("Start to QC SNPs for combined file of gwas1 and gwas2")
  # matched SNP but unmatched allele
  raw_n = nrow(gwas12)
  gwas12 = gwas12[!is.na(b_2), ]
  del_n = raw_n - nrow(gwas12)
  if(del_n > 0) {
    log_info(sprintf("There were %d SNPs have unmatched alleles between gwas1 and gwas2, will delete these SNPs!", del_n))
  } else {log_info("All SNPs in merge file were matched alleles")}
  
  # freq difference 
  raw_n = nrow(gwas12)
  gwas12 = gwas12[abs(freq_1 - freq_2) <= diffreq, ]
  del_n = raw_n - nrow(gwas12)
  if(del_n > 0){
    log_info(sprintf("There were %d SNPs have different AF > %.2f, will delete these SNPs!", del_n, diffreq))
  } else {log_info(sprintf("All SNPs different freq in merge file were < %.2f", diffreq))}

  # raw_n = nrow(gwas12)
  # gwas12 = gwas12[b_1 * b_2 > 0, ]
  # del_n = raw_n - nrow(gwas12)
  # if(del_n > 0){
  #   log_info(sprintf("There were %d SNPs have unmatched effect direction, will delete these SNPs!", del_n))
  # } else {log_info("All SNPs effect size in merge file were matched direction")}
  
  log_info("Start to demeta SNP effect size and se")
  gwas12[, denom := 1/se_1^2 - 1/se_2^2]
  anom = subset(gwas12, gwas12$denom < 0)
  log_info(sprintf("There were %d SNPs that se smaller in gwas2!", nrow(anom)))

  se_mean = mean(gwas12$se_1, na.rm=TRUE)
  gwas12[, se_demeta := se_mean]
  gwas12[denom > 0, se_demeta := sqrt(1/denom)]
  gwas12[, denom := NULL]
  
  # gwas12[, se_demeta := fifelse((1/se_1^2 - 1/se_2^2) > 0, sqrt(1/(1/se_1^2 - 1/se_2^2)), se_mean)]
  gwas12[, b_demeta := b_1-se_demeta^2 * (b_2-b_1)/se_2^2]
  
  raw_n = nrow(gwas12)
  gwas12 = gwas12[b_demeta != 0, ]
  del_n = raw_n - nrow(gwas12)
  if (del_n > 0) {
    log_info(sprintf("There were %d SNPs has been deleted casued by the demeta value is zero!", del_n))
  }

  # het test 
  gwas12[, z_het := (b_demeta - b_2) / sqrt(se_demeta^2 + se_2^2)]
  gwas12[, p_het := 2*pnorm(-abs(z_het))]
  
  log_info("Start to create effect N of demeta")  
  gwas12[, N := N_1 - N_2]
  raw_n = nrow(gwas12)
  gwas12 = gwas12[N > 0, ]
  del_n = raw_n - nrow(gwas12)
  if (del_n > 0) {
    log_info(sprintf("There were %d SNPs has been deleted casued by the demeta N < 0!", del_n))
  }

  outfile = gwas12[, c("SNP","A1_1","A2_1","freq_1","b_demeta","se_demeta","p_1","N", "z_het", "p_het")]
  names(outfile) = c("SNP","A1","A2","freq","b","se","p","N", "Het", "p_het")
  
  log_info("Start to calculate P value based on beta and se")
  outfile$p <- 2 * (1 - pnorm(abs(outfile$b / outfile$se)))
  
  log_info("Demeta analysis is completed!") 
  return(outfile)
}

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/GWAS")
library(data.table)
meta = fread("HP.sex.combined.GWAS.tsv.gz", select=c("rs_id","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","n"))
sub_c = fread("2247_1.v1.0.fastGWA.gz", select=c("SNP","A1","A2","AF1","BETA","SE","P","N"))

gwas1 = gwasQC(meta, MAF=0.01)
gwas2 = gwasQC(sub_c, MAF=0.01)
gwas2[, `:=` (b  = round(b, 4), se = round(se, 4))]

res = de_meta(gwas1, gwas2)
res_het = res[p_het > 0.001, ]
res_het1 = merge(res_het, gwas1, by="SNP", suffixes = c("_de", "_m"))
res_het2 = res_het1[!(abs(b_de) > 1 | (p_de < 5e-8 & p_m > 5e-8)), ][order(-abs(b_de))]

info = fread("HP.sex.combined.GWAS.tsv.gz", select=c("chromosome", "base_pair_location", "rs_id"))
res2= merge(res_het2, info, by.x="SNP", by.y="rs_id")[, c("chromosome", "base_pair_location", "SNP", "A1_de", "A2_de", "freq_de", "b_de", "se_de", "p_de", "N_de")]
names(res2) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")

intercept = 1.2087
res2[, z_int :=  beta/SE / sqrt(intercept)]
res2[, se_int := beta / z_int]
res2[, p_int := 2*pnorm(-abs(z_int))]

res2 = res2[, .(CHR, POS, SNP, A1, A2, freq, beta, se_int, p_int, N)]
names(res2) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")
fwrite(res2, file="Demeta_HP.sex.combined.rawN_Het0001_fltP_intercept.GWAS.gz", sep="\t")


#==========================================================#
#       *^  METAL for meta analysis  ^* 
#==========================================================#
#============ The EAS meta BBJ and MVP ============#
# METAL Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
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
GENOMICCONTROL OFF

ADDFILTER freq > 0.01
ADDFILTER freq < 0.99
TRACKPOSITIONS ON 
CHROMOSOME CHR
POSITION POS
# additional options
EFFECT_PRINT_PRECISION 12
STDERR_PRINT_PRECISION 12
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/BBJ/hearing-loss/BBJ_EAS_reformat.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_EAS_reformat_hg19.gz

OUTFILE EAS_MVP_BBJ .tbl
ANALYZE HETEROGENEITY 
QUIT
#=========================
# run metal
metal="/public/home/shilulu/software/METAL/build/bin/metal"
qsubshcom "$metal EAS_MVP_BBJ_metal.conf" 1 100G METAL 1:00:00 ""
#=========================

#============ The EUR meta Trpchevska et al., MVP and demeta of De Angelis et al. ============#
# METAL Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
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
GENOMICCONTROL OFF

ADDFILTER freq > 0.01
ADDFILTER freq < 0.99
TRACKPOSITIONS ON 
CHROMOSOME CHR
POSITION POS
# additional options
EFFECT_PRINT_PRECISION 12
STDERR_PRINT_PRECISION 12
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/UKB/hearing-loss/AJHG_EUR_reformat.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_EUR_reformat_hg19.gz
PROCESS /public/home/shilulu/Wulab_project/ARHL/NC_sup_test/00.GWAS/Demeta_HP.sex.combined.rawN_Het0001_fltP_intercept.GWAS.gz

OUTFILE EUR_MVP_Trpchevska_De-Angelis .tbl
ANALYZE HETEROGENEITY 
QUIT
#=========================
# run metal
metal="/public/home/shilulu/software/METAL/build/bin/metal"
qsubshcom "$metal EUR_MVP_Trpchevska_De-Angelis_metal.conf" 1 100G METAL 1:00:00 ""
#=========================


#======================
# reformat and filter METAL result 
# the EAS and EUR, which filter the SNP only occured in one study
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
#=============================#
#  meta for EAS EUR AFR AMR   #
# the mean result to analysis #
#=============================#
# METAL Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
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
GENOMICCONTROL OFF

ADDFILTER freq > 0.01
ADDFILTER freq < 0.99
TRACKPOSITIONS ON 
CHROMOSOME CHR
POSITION POS
# additional options
EFFECT_PRINT_PRECISION 12
STDERR_PRINT_PRECISION 12
CUSTOMVARIABLE TotalSampleSize  
LABEL TotalSampleSize as N

PROCESS /public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/EUR_MVP_Trpchevska_De-Angelis_METALfilter.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_AFR_reformat_hg19.gz
PROCESS /public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss/MVP_AMR_reformat_hg19.gz
PROCESS /public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/EAS_MVP_BBJ_METALfilter.gz

OUTFILE All_MVP_Trpchevska_De-Angelis_BBJ .tbl
ANALYZE HETEROGENEITY 
QUIT 

#=========================
# run metal
metal="/public/home/shilulu/software/METAL/build/bin/metal"
qsubshcom "$metal All_MVP_Trpchevska_De-Angelis_BBJ_metal.conf" 1 100G METAL 90:00:00 ""

#==================================
# reformat and filter METAL result 
# and add chr pos for the result
#  rfm_flt_addchrpos_metal.r
#==================================
library(dplyr)
library(stringr)
library(data.table)

setwd("/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/01.meta")
df <- fread( "All_MVP_Trpchevska_De-Angelis_BBJ1.tbl")[, c("Chromosome", "Position", "MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize")]
names(df) = c("CHR", "POS", "SNP", "A1", "A2", "freq", "beta", "SE", "p", "HetISq", "HetChiSq", "HetDf", "HetPVal", "N")

# df <- df %>%
#     filter(HetDf > 0) %>%
#     mutate(A1 = str_to_upper(A1), A2 = str_to_upper(A2))
# fwrite(df, "h2/All_MVP_Trpchevska_De-Angelis_BBJ.txt.gz", sep="\t")
# filter the SNP which only occured at one study
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
##########
qsubshcom "Rscript rfm_flt_addchrpos_metal.r" 1 100G METAL 1:00:00 ""
#############
