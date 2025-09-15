#### 2. forest plot for GSMR
setwd("F:/Github/PHD_job/2_project_hearing loss/NC_revision/badgers/gsmr")

library(ggplot2)
library(data.table)
library(patchwork)

data = fread("ARHL_gsmr.gsmr")
data$FDR <- p.adjust(data$p, method='BH')
final = data[FDR < 0.05, ][, 1:7]
fwrite(data, "ARHL_gsmr_noNA.FDR.gsmr", sep="\t")
fwrite(final, "ARHL_gsmr_noNA.FDR005.gsmr", sep="\t")

final$lower=final$bxy-1.96*final$se
final$upper=final$bxy+1.96*final$se

data1 = final[Exposure != "ARHL",]
data1[data1$Exposure=="Tiredness","Exposure"]="Frequency of tiredness / lethargy in last 2 weeks"
data1[data1$Exposure=="Neuroticism","Exposure"]="Neuroticism score"
data1[data1$Exposure=="LongStandIllness","Exposure"]="Long-standing illness"
data1[data1$Exposure=="Overall_health_rating","Exposure"]="Overall health rating"
data1[data1$Exposure=="derpess_mood","Exposure"]="Frequency of depressed mood in last 2 weeks"
data1[data1$Exposure=="Past_tobacco_smoking","Exposure"]="Past tobacco smoking"

data1$Exposure = factor(data1$Exposure, levels=rev(data1$Exposure))
p1 = ggplot(data=data1, aes(y=Exposure, x=bxy)) +
  geom_point( shape=19, size=3, color="#1B9E77") +
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

ggsave("gsmr_plot.png", p1, dpi=600, width=5.5, height=3)
ggsave("gsmr_plot.pdf", p1, width=5.5, height=3)