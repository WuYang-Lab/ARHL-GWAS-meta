# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# replication using finngen Sensorineural hearing loss
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta/replication")
library(data.table)
library(dplyr)

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

lead_snp = fread("/public/home/shilulu/project_hearing-loss/new_run/all_meta/ld_clump/ARHL_MVP_AJHG_BBJ_reformatMETAL.clumped.lead.snp")
gwas_meta = fread("/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL")
gwas_finngen = fread("finngen_R11_H8_HL_SEN_NAS.gz")[, c("rsids", "ref", "alt", "beta", "sebeta", "af_alt", "pval")]
# alt is effect allele
names(gwas_finngen) = c("SNP", "ref", "alt", "b", "se", "freq", "pval")
gwas_finngen = filter(gwas_finngen, !is.na(SNP) & SNP != "")
gwas_finngen$N = 39620 + 397865
gwas_finngen_rescale = rescale_beta_se(gwas_finngen, 0.1, 39620, 397865)
names(gwas_meta) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
gwas_meta_rescale = rescale_beta_se(gwas_meta, 0.1, 438594, 1035810)

lead_meta = left_join(lead_snp, gwas_meta_rescale, by="SNP") %>%
  filter(!is.na(b))
lead_finngen = left_join(lead_snp, gwas_finngen_rescale, by="SNP") %>%
  filter(!is.na(b))
fwrite(lead_finngen, file="lead_snp_in_finngen.txt", sep="\t")


novel_snp = fread("novel_snp_ARHL.txt")
meta_finngen_snp = left_join(meta_finngen, novel_snp, by="SNP") %>% 
  mutate(cate = if_else(is.na(CHR | BP), "Known", "Previously unknown")) %>%
  select(-c(CHR, BP))

cor_file = cor.test(meta_finngen_snp$b1, meta_finngen_snp$b2, method="pearson")
p_value = cor_file$p.value
cor_value = cor_file[["estimate"]][["cor"]]
meta_finngen_snp$y_upper = cor_file$conf.int[2] * meta_finngen_snp$b1
meta_finngen_snp$y_lower = cor_file$conf.int[1] * meta_finngen_snp$b1
meta_finngen_snp$x_max = meta_finngen_snp$b1 + meta_finngen_snp$se1 * 1.96
meta_finngen_snp$x_min = meta_finngen_snp$b1 - meta_finngen_snp$se1 * 1.96
meta_finngen_snp$y_max = meta_finngen_snp$b2 + meta_finngen_snp$se2 * 1.96
meta_finngen_snp$y_min = meta_finngen_snp$b2 - meta_finngen_snp$se2 * 1.96

p_rep = ggplot(meta_finngen_snp, aes(x=b1, y=b2, color=cate)) +
  geom_errorbar(aes(xmin=x_min, xmax=x_max), width=0, color="grey") +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0, color="grey")+
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.8) +
  geom_vline(xintercept=0, linetype="dashed", color="black", size=0.8) +
  geom_abline(slope=cor_value, intercept=0, color="red", size=0.8)+
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = "grey", color=NA, alpha = 0.5) +
  geom_point() +
  scale_color_manual(values=c("Previously unknown" = "#E3191C", "Known" = "#92469E"))+
  labs(x="ARHL Meta", y="HL Finngen", 
       title=paste("Correlation = ", round(cor_value, 3), ", ", "P-value = ", sprintf("%.2e", p_value))) +
  theme_classic() +
  theme(legend.position=c(0.8,0.15),
        legend.title = element_blank())
p_rep

ggsave("replication_pearson_cor.png", p_rep, height=5, width=5, dpi=500)

save(meta_finngen, lead_meta, lead_finngen, file="Finngen_replication.RData")
###########################################################################