#- - - - - demeta UKB effect 
rescale_b_se <- function(gwas, case, control){
  names(gwas) = c("SNP","A1","A2","freq","b","se","p","N")
  n_eff       = as.integer(4*case*control/(case + control))
  p           = gwas$freq
  z           = gwas$b / gwas$se
  se          = 1/sqrt(2*p*(1-p)*(n_eff+z^2))
  b           = z*se
  gwas$b      = b
  gwas$se     = se
  gwas$N      = n_eff
  
  return(gwas)
}


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

setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/Meta")
library(data.table)
library(ggthemes)
meta = fread("All_MVP_Trpchevska_De-Angelis_BBJ_filter.gz")
names(meta) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
re_meta = rescale_b_se(meta, 456613, 1053834)

sub_c = fread("2247_1.v1.0.fastGWA.gz", select=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
names(sub_c) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
re_sub_c = rescale_b_se(sub_c, 114318, 323449) 

gwas1 = gwasQC(re_meta, MAF=0.01)
gwas2 = gwasQC(re_sub_c, MAF=0.01)

res = de_meta(gwas1, gwas2)
res_het = res[p_het > 0.001, ]
res_het1 = merge(res, meta, by="SNP", suffixes = c("_de", "_m"))
res_het2 = res_het1[!(abs(b_de) > 1 | (p_de < 5e-8 & p_m > 5e-8)), ][order(-abs(b_de))]

#- hm3
hm3 <- fread("listHM3.txt", header = FALSE, col.names = "SNP")
res_hm3 = merge(res_het1, hm3, by="SNP")
res2 = res_hm3[, .(SNP, A1_de, A2_de, freq_de, b_de, se_de, p_de, N_de)]
names(res2) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

fwrite(res2, "Demeta_raw_hm3_Meta_rescale.txt", sep="\t")


p <- ggplot(res_hm3, aes(y = b_de / se_de, x = b_m / se_m)) +
  geom_point(color = "#69b3a2") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(y=expression(Zscore[DE_META]), x = expression(Zscore[META])) +
  theme_base() 

ggsave("badgers_demetaVmeta.png", p, width=4, height=4, dpi=600)