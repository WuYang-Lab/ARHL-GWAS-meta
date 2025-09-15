#- - - - - - - multi MR
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript
source("/public/home/shilulu/script/differentMR/gsmr_plot_function.R")
source("/public/home/shilulu/script/differentMR/GSMR_effect_plot.R")
source("/public/home/shilulu/script/differentMR/comfun_file.R")
source("/public/home/shilulu/script/differentMR/MR.R")

library(survey)
library(NCmisc)
library(MASS)

# ex_out must have 2 columns which contain the Exposure and Outcome (the significant results, you know!)
setwd("/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR")
ex_out <- read.table("ARHL_gsmr_noNA.FDR005.gsmr", header=TRUE, sep="\t")
gsmr_data1 = read_gsmr_data("ARHL_gsmr_noheidi.eff_plot.gz")

res=c()
for (k in 1:nrow(ex_out)) {
  ex <- ex_out$Exposure[k]
  out <- ex_out$Outcome[k]

  resbuf = gsmr_snp_effect(gsmr_data1, ex, out)
  bzx=resbuf$bzx;bzx_se=resbuf$bzx_se;bzx_pval=resbuf$bzx_pval
  bzy=resbuf$bzy;bzy_se=resbuf$bzy_se;bzy_pval=resbuf$bzy_pval

  ## Main MR function    
  # for(method in c("GSMR_heidi_v1","GSMR_heidi_v3_stepwise", "GSMR_noheidi", "IVW", "Median", "Simple_Median", "Mode", "Egger", "Robust", "Lasso", "RAPS", "PRESSO", "MRMix", "Con-Mix")){
  for(method in c("GSMR_heidi_v1", "IVW", "Simple_Median", "Robust", "RAPS", "PRESSO")){
      res_mr = all_mr(bzx=resbuf$bzx, bzx_se=resbuf$bzx_se, bzx_pval=resbuf$bzx_pval,
                      bzy=resbuf$bzy, bzy_se=resbuf$bzy_se, bzy_pval=resbuf$bzy_pval, method=method)
      res = rbind(res,data.frame(exposure=ex, outcome=out, method, p=res_mr$pval, bxy=res_mr$bxy,nr_pleio_rm=length(res_mr$pleio)))
    }
  print(paste(ex,"to",out,"completed",sep = " "))
}

res$z=p.to.Z(res$p)
res$se=res$bxy/res$z
res$z=res$z*sign(res$se)
res$se=res$se*sign(res$se)
res$association=paste(res$exposure,"_to_",res$outcome,sep="")
res$lower=res$bxy-1.96*res$se
res$upper=res$bxy+1.96*res$se

setwd("/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR")
write.table(res, file="ARHL_MR_comparison_results.txt", row=F, col=T, quote=F)

#- - - - - - - - - plot 
setwd("F:/Github/PHD_job/2_project_hearing loss/NC_revision/badgers/gsmr")

library(data.table)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

resdt=read.table("ARHL_MR_comparison_results.txt", header=T)
#  use the GSMR result you fig to pheatmap
gsmr_gene = read.table("ARHL_gsmr_noNA.FDR005.gsmr", header=T)

res <- resdt |> 
  left_join(gsmr_gene, by = c("exposure" = "Exposure", "outcome" = "Outcome")) |> 
  mutate(bxy = if_else(method == "GSMR_heidi_v1", bxy.y, bxy.x),
         se = if_else(method == "GSMR_heidi_v1", se.y, se.x),
         p = if_else(method == "GSMR_heidi_v1", p.y, p.x)) |> 
  select(-bxy.y, -se.y, -p.y,
         -bxy.x, -se.x, -p.x)

res$lower=res$bxy-1.96*res$se
res$upper=res$bxy+1.96*res$se

res[res$method=="GSMR_heidi_v1","method"]="GSMR"
res[res$method=="Simple_Median","method"]="Simple Median"
res[res$exposure=="Tiredness","exposure"]="Frequency of tiredness / lethargy in last 2 weeks"
res[res$exposure=="Neuroticism","exposure"]="Neuroticism score"
res[res$exposure=="LongStandIllness","exposure"]="Long-standing illness"
res[res$exposure=="Overall_health_rating","exposure"]="Overall health rating"
res[res$exposure=="derpess_mood","exposure"]="Frequency of depressed mood in last 2 weeks"
res[res$exposure=="Past_tobacco_smoking","exposure"]="Past tobacco smoking"
res[res$exposure=="Waist_circumference","exposure"]="Waist circumference"


res$method <- factor(res$method, levels = levels(factor(res$method)))
# res$association = factor(res$association, levels=levels(res$association))

colors <- rev(brewer.pal(8, "Set2"))
library(patchwork)
p <- ggplot(data=res, aes(y=method,x=bxy,color=method))+
  geom_point(aes(color=method),shape=19, size=3)+
  geom_errorbar(aes(xmin=lower,xmax=upper),width=0.3,size=1)+
  geom_vline(xintercept=0,linetype="dashed") +
  
  scale_color_manual(values=colors) +
  theme_few() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, angle=60, hjust=1, vjust=1),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=12))

plots <- lapply(split(res, res$exposure), function(df) {
  p %+% df + labs(x="", y="", title = unique(df$exposure))
})

plots[[1]] <- plots[[1]] + 
  theme(axis.text.y = element_text(size=12))
plots[[5]] <- plots[[5]] + 
  theme(axis.text.y = element_text(size=12))

p2 <- wrap_plots(plots, nrow = 2)
ggsave("Multi_method_MR.png", p2, width=9, height=6.5, dpi=600)
ggsave("Multi_method_MR.pdf", p2, width=9, height=6.5)
