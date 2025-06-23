#==============================#
#   plot Genetic architecture  #
#==============================#
# plot MAF vs effect size
# select clump lead SNP to plot
setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/plot_Genetic_architecture")

library(dplyr)
library(data.table)
# plot MAF and effect size, recommand use clumped SNP file 
data = read.table("plot_mat_vs_beta.txt", sep="\t", header=TRUE)
data$MAF = if_else(data$freq < 0.5, data$freq, 1 - data$freq)
data$b = if_else(data$freq < 0.5, data$beta, -data$beta)
data$p1  = -log10(data$p)

library(ggplot2)
library(viridis)
p_maf = ggplot(data, aes(x=MAF, y=b, color=p1)) +
  geom_point(alpha=0.7, size=2) +
  theme_classic() +
  scale_color_viridis() +
  labs(x="Minor Allele Frequency (MAF)", 
       y="Effect Size (BETA)",
       color=expression(-log[10](P-value)))+
  theme(legend.position = c(0.9, 0.8),
        legend.title.position = "left",
        legend.title = element_text(angle=90))+
  theme(axis.text.x=element_text(size=12, hjust=1, vjust=1),
        axis.text.y=element_text(size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))
p_maf

setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/plot_Genetic_architecture")
plt_list = data.table()
for (chr in 1:22){
  res = fread(paste0("bayesS_annot_chr/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesS_chr", chr, ".parRes"))
  res$chr = chr
  plt_list = rbind(plt_list, res[V1=="Pi" | V1=="S" | V1=="hsq", ])
}
colnames(plt_list)[1] = c("cate")
# fwrite(plt_list, "gctb_result.txt", sep="\t")

plt_list$y_max = plt_list$Mean + plt_list$SD
plt_list$y_min = plt_list$Mean - plt_list$SD

chr_len=fread("chr_length.txt")
plt_dt = merge(plt_list, chr_len, by="chr")

h = sum(subset(plt_list, cate=="hsq")$Mean)
pim = mean(subset(plt_list, cate=="Pi")$Mean)
dt = subset(plt_list, cate=="S")
dt[, weight := 1 / (SD^2)]
S_combined <- sum(dt$Mean * dt$weight) / sum(dt$weight)
SE_combined <- sqrt(1 / sum(dt$weight))
> S_combined
[1] -0.5297184
> SE_combined
[1] 0.01920526

# pearson correaltion for hsq
h_dt = subset(plt_dt, cate=="hsq")
pearson_cor = cor(h_dt$Mean, h_dt$bp, method="pearson")
pearson_cor
# [1] 0.9134548

s = ggplot(subset(plt_list, cate=="S"), aes(x=chr, y=Mean)) +
  geom_bar(stat="identity", position="dodge", fill="#B03060") +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  scale_x_continuous(breaks = 1:22) +
  labs(x="CHR", y="S")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12,angle=90, hjust=1, vjust=1),
        axis.text.y=element_text(size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))

hsq = ggplot(subset(plt_dt, cate=="hsq"), aes(x=bp, y=Mean)) +
  geom_smooth(aes(x=bp, y=Mean), method="lm", level=0.95)+
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  geom_point(aes(y = Mean), size = 0.5) +
  geom_text(aes(y = Mean, label = chr), hjust = -0.2, vjust = -0.5, size = 3) + 
  labs(x="Chromosome length (Mb)", y="Heritability")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, hjust=1, vjust=1),
        axis.text.y=element_text(size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))

pi = ggplot(subset(plt_list, cate=="Pi"), aes(x=chr, y=Mean)) +
  geom_bar(stat="identity", position="dodge", fill="#018B00") +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  scale_x_continuous(breaks = 1:22) +
  scale_y_continuous(expand=c(0, 0)) + 
  labs(x="CHR", y="Pi")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=1),
        axis.text.y=element_text(size=12),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))



library(patchwork)
p_all = (p_maf | s) / (hsq | pi) + 
  plot_annotation(tag_levels = 'a', tag_prefix='', tag_suffix='')

p_all

setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/plot_Genetic_architecture")
ggsave("plot_genetic_architecture.png", p_all, height=7, width=9, dpi=500)
ggsave("plot_genetic_architecture.pdf", p_all, height=7, width=9)
