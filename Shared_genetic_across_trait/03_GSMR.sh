#--------------------------- GSMR ----//
#
# GSMR for some traits select from BADGERS
# 
#--------------------------- GSMR ----//

###################################################
#------------- QC for GWAS summary ----//
sums_col_treate <- function(gwas, SNP, A1, A2, freq, P,
                            CHR = NULL, POS = NULL,
                            BETA = NULL, SE = NULL,
                            OR = NULL, Z = NULL, N=NULL) {
    gwas_col = c(
        CHR   =   CHR,
        POS   =   POS, 
        SNP   =   SNP,
        A1    =    A1, 
        A2    =    A2,
        freq  =  freq,
        b     =  BETA,
        se    =    SE,
        OR    =    OR,
        Z     =     Z,
        p     =     P,
        N     =     N)
    gwas_col = gwas_col[!sapply(gwas_col, is.null)]

    missing_required = setdiff(c("SNP", "A1", "A2", "freq", "p"), names(gwas_col))
    if (length(missing_required) > 0) {
        stop("Missing required columns: ", paste(missing_required, collapse = ", "))
    }
    has_b_se = all(c("b", "se") %in% names(gwas_col))
    has_or   = "OR" %in% names(gwas_col)
    has_z    = "Z" %in% names(gwas_col)
    if (!(has_b_se || has_or || has_z)) {
        stop("Must provide at least one of the following: (BETA + SE), OR, or Z.")
    }

    gwas = gwas[, unname(gwas_col), with=FALSE]
    setnames(gwas, old = unname(gwas_col), new = names(gwas_col))
    
    valid_p <- !is.na(gwas$p) & gwas$p >= 0 & gwas$p <= 1
    if (has_b_se){    
        gwas = gwas[valid_p & !is.na(b) & !is.na(se) & b != 0, ]
    } else if (has_or){
        gwas = gwas[valid_p & !is.na(OR), ]
    } else if (has_z){
        gwas = gwas[valid_p & !is.na(Z), ]
    }
    
    return(gwas)
}

library(data.table)
gwas = fread("Ever_smoke.fastGWA.gz")
clean = sums_col_treate(gwas, "variant_id", "effect_allele", "other_allele", "effect_allele_frequency", "p_value", BETA="beta", SE="standard_error" N="N")
clean = sums_col_treate(gwas, "SNP", "A1", "A2", "AF1", "P", BETA="BETA", SE="SE", N="N")
clean = clean[, .(SNP, A1, A2, freq, b, se, p, N)]
fwrite(clean, "../cleanGWAS/clean_ever_smoke.fastGWA.gz", sep="\t")

#############################################
dir="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS"
echo "Loneliness $dir/clean_Loneliness_GCST006924.txt.gz
Insomnia $dir/clean_Sleeplessness-insomnia.fastGWA.gz
LongStandIllness $dir/clean_Long-standing_illness.fastGWA.gz
Overall_health_rating $dir/clean_Overall_health_rating.fastGWA.gz
Snoring $dir/clean_Snoring.fastGWA.gz
Neuroticism $dir/clean_Neuroticism_score.fastGWA.gz
Tiredness $dir/clean_tiredness_lethargy2080.fastGWA.gz
Leg_fat_percentage $dir/clean_leg_fat_percentage_right.fastGWA.gz
Noisy_workplace $dir/clean_Noisy_workplace.fastGWA.gz
Taking_other_prescription_medications $dir/clean_Taking_other_prescription_medications.fastGWA.gz
derpess_mood $dir/clean_depressed_mood.fastGWA.gz
Ever_smoke $dir/clean_ever_smoke.fastGWA.gz" > exposure_new.txt

echo "ARHL /public/home/shilulu/Wulab_project/ARHL/14May2025/Demeta_raw_hm3_Meta_rescale_keeponly_SNP_14May2025.txt" > outcome.txt

#=============================================#
#! /bin/bash 
gcta="/public/home/wuyang/bin/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
ref="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
exposure="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS/exposure_new.txt"
outcome="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS/outcome.txt"

cmd="${gcta} --bfile ${ref} \
--gsmr-file ${exposure} ${outcome} \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--diff-freq 0.5 \
--clump-r2 0.05 \
--gsmr2-beta \
--gsmr-snp-min 10 \
--heidi-thresh 0.01 \
--effect-plot \
--out result/ARHL_gsmr_new \
--thread-num 20"
qsubshcom "$cmd" 20 100G gsmr 90:00:00 ""


#========================================================#
### annote the IDP region, also delete row contained "nan"
#========================================================#
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript
library(data.table)
library(dplyr)

gsmr <- fread("ARHL_gsmr_new.gsmr")
gsmr[gsmr == "nan"] <- NA

# p < 0.05 / row number / 2 for 2 dirction gsmr
# remove row contained "nan" (all)
final <- gsmr %>%
  filter(p < 0.05 / (nrow(gsmr))) %>%
  filter(complete.cases(.))

fwrite(final, file="ARHL_gsmr_new_noNA_anno.gsmr", sep="\t", na="nan", quote=FALSE)

#### 2. forest plot
setwd("F:/Github/PHD_job/2_project_hearing loss/process/badgers/gsmr")

library(ggplot2)
library(data.table)
library(patchwork)
data = fread("ARHL_gsmr_new.gsmr")

data$lower=data$bxy-1.96*data$se
data$upper=data$bxy+1.96*data$se


data$FDR <- p.adjust(data$p, method='BH')
final = data[FDR < 0.05, ]

final$fill = ifelse(final$p < 0.05 / 24, "type1", "type2")

final$Exposure = factor(final$Exposure, levels = rev(unique(final$Exposure)))
final$Outcome = factor(final$Outcome, levels = rev(unique(final$Outcome)))

data1 = final[Exposure != "ARHL",]
data2 = final[Exposure == "ARHL",]


p1 = ggplot(data=data1, aes(y=Exposure, x=bxy)) +
  geom_point(aes(fill = fill), shape=21, size=3, color="#1B9E77") +
  scale_fill_manual(values = c(type1 = "#1B9E77", type2 = "white")) + 
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.2, linewidth=0.8, color="#1B9E77")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.background = element_blank(),
        axis.title = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "none")
p1

p2 = ggplot(data=data2, aes(y=Outcome, x=bxy)) +
  geom_point(aes(fill = fill), shape=21, size=3, color="#D95F02") +
  scale_fill_manual(values = c(type1 = "#D95F02", type2 = "white")) + 
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.2, linewidth=0.8, color="#D95F02")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + 
  theme(axis.text.x = element_text( size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.background = element_blank(),
        axis.title = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "none")
p2


p = p1 + p2 + plot_layout(widths = c(2, 2))
p

ggsave("gsmr_plot.png", p, dpi=600, width=9, height=3)
ggsave("gsmr_plot.pdf", p, width=9, height=3)