#===================================================#
#       Deming regression | pearson correlation     #
#         GWAS significant clumped
#===================================================#
# # # rescale beta and se
rescale_beta_se <- function(df, K, case, control){
  K = K # pop prevalence
  i = dnorm(qnorm(1-K))/K
  v = case / (case + control) # case / total, sample prevalence 
  n = df$N 
  n_eff = i^2*v*(1-v)/(1-K)^2 * n
  freq = df$freq 
  z = df$b / df$se
  # new beta and se
  se = 1/sqrt(2*freq*(1-freq)*(n_eff+z^2))
  b = se*z

  df$b = b 
  df$se = se 
  df$N = n_eff
  return(df)
}

setwd("/public/home/shilulu/project_hearing-loss/new_run")
library(data.table)
df_eur <- fread("ARHL_EUR_MVP_AJHG_reformatMETAL.gz")[, c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")]
setnames(df_eur, c("beta", "SE"), c("b", "se"))
EUR_df = rescale_beta_se(df_eur, 0.1, 377577, 751989)

df_eas <- fread("ARHL_EAS_MVP_BBJ_reformatMETAL.gz")[, c("SNP", "A1", "A2", "freq", "beta", "SE", "p", "N")]
setnames(df_eas, c("beta", "SE"), c("b", "se"))
EAS_df = rescale_beta_se(df_eas, 0.1, 6071, 178628)

setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/ARHL_meta")
fwrite(EUR_df, file="EUR_hearing_loss_rescale_betase.gz", sep="\t", row.names=FALSE, quote=FALSE)
fwrite(EAS_df, file="EAS_hearing_loss_rescale_betase.gz", sep="\t", row.names=FALSE, quote=FALSE)
#############################################
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss")
library(dplyr)
library(data.table)
df <- fread("MVP_R4.1000G_AGR.Phe_389.AMR.GIA.dbGaP.txt.gz")
# rename columns for mrmega in format
old_col <- c("SNP_ID", "chrom", "pos", "ref", "alt", "af", "num_samples", "or", "pval")
new_col <- c("SNP", "CHR", "BP", "A2", "A1", "freq", "N", "OR", "P")
setnames(df, old_col, new_col)

# OR to BETA
df$b <- log(df$OR)
# if beta=0 (represent OR=1 and P=1), then se=mean(allse)
se_mean <- df %>%
  filter(b!=0) %>%
  mutate(se = sqrt(b^2 / qchisq(P, 1, lower.tail=FALSE))) %>%
  summarise(mean_se = mean(se, na.rm = TRUE)) %>%
  pull(mean_se)

df$se <- ifelse(df$b==0, se_mean, sqrt(df$b^2/qchisq(df$P, 1, lower.tail=F)))

## rescale beta and se
df_amr = rescale_beta_se(df, 0.1, 23888, 28460)
# save the file
out_col <- df_amr[, .(SNP, CHR, BP, A1, A2, freq, b, se, P, N)]
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/ARHL_meta")
fwrite(out_col, file="MVP_AMR_hearing_loss_rescale_betase.txt", sep="\t", row.names=FALSE)
#############################################

#############################################
### data rescale and reformat for lead SNP based correlation 
library(dplyr)
library(data.table)

extract_lead <- function(gwas){
  lead_meta = fread("/public/home/shilulu/project_hearing-loss/new_run/all_meta/ld_clump/ARHL_MVP_AJHG_BBJ_reformatMETAL.clumped.lead.snp")
  pop = fread(gwas)[, c("SNP", "A1", "A2", "b", "se")]
  # extract lead SNP beta from ancestry respectively
  lead_pop = left_join(lead_meta, pop, by="SNP") %>%
    filter(!is.na(b)) %>%
    select(SNP, A1, A2, b, se)
  return(lead_pop)
}

fix_effect <- function(pop1, pop2){
  # fixed the beta effect SNP
  df = inner_join(pop1, pop2, by="SNP") %>%
    mutate(b2 = case_when(A1.x == A1.y & A2.x == A2.y ~ b.y,
                          A1.x == A2.y & A2.x == A1.y ~ -b.y)) %>% 
    rename(b1 = b.x, se1 = se.x, se2 = se.y) %>% 
    select(b1, se1, b2, se2)
  return(df)
}

run_pearson <- function(df){
  correlation_pearson = cor.test(df$b1, df$b2, method="pearson")
  return(correlation_pearson)
}

plot_cor_pearson <- function(df, cor_file, labx, laby){
  p_value <- cor_file$p.value
  cor_value <- cor_file[["estimate"]][["cor"]]
  df$y_upper <- cor_file$conf.int[2] * df$b1
  df$y_lower <- cor_file$conf.int[1] * df$b1
  df$x_max = df$b1 + df$se1 * 1.96
  df$x_min = df$b1 - df$se1 * 1.96
  df$y_max = df$b2 + df$se2 * 1.96
  df$y_min = df$b2 - df$se2 * 1.96
  
  library(ggplot2)
  p=ggplot(df, aes(x=b1, y=b2)) +
    geom_errorbar(aes(xmin=x_min, xmax=x_max), width=0, color="grey") +
    geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0, color="grey")+
    geom_hline(yintercept=0, linetype="dashed", color="black", size=0.8) +
    geom_vline(xintercept=0, linetype="dashed", color="black", size=0.8) +
    geom_abline(slope=cor_value, intercept=0, color="red", size=0.8)+
    geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "grey", alpha = 0.5) +
    geom_point(color="#1F78B4") +
    labs(x=labx, y=laby, 
         title=paste("Correlation = ", round(cor_value, 3), ", ", "P-value = ", sprintf("%.2e", p_value))) +
    theme_classic()
  return(p)
}

################################################################################################
setwd("/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/ARHL_meta")
eur_lead = extract_lead("EUR_hearing_loss_rescale_betase.gz")
eas_lead = extract_lead("EAS_hearing_loss_rescale_betase.gz")
amr_lead = extract_lead("MVP_AMR_hearing_loss_rescale_betase.txt")
afr_lead = extract_lead("MVP_AFR_hearing_loss_rescale_betase.txt.gz")
save(eur_lead, eas_lead, amr_lead, afr_lead, file="lead_SNP_in_ancesatry.RData")
save.image("for_Deming_Pearson.RData")
################################################################################################

eur_eas_cor = run_pearson(fix_effect(eur_lead, eas_lead))
eur_afr_cor = run_pearson(fix_effect(eur_lead, afr_lead))
eur_amr_cor = run_pearson(fix_effect(eur_lead, amr_lead))


df1 = fix_effect(eur_lead, eas_lead)
df1_cor = run_pearson(df1)
p1 = plot_cor_pearson(df1, df1_cor, "EUR", "EAS")
p1

df2 = fix_effect(eur_lead, afr_lead)
df2_cor = run_pearson(df2)
p2 = plot_cor_pearson(df2, df2_cor, "EUR", "AFR")
p2

df3 = fix_effect(eur_lead, amr_lead)
df3_cor = run_pearson(df3)
p3 = plot_cor_pearson(df3, df3_cor, "EUR", "AMR")
p3


library(patchwork)
p_all <- p1 | p2 | p3
p_all
ggsave("pearson_cor.png", p_all, height=4, width=12, dpi=500)