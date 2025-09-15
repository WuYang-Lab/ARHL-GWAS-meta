# # # process result file and plot
library(data.table)
library(ggplot2)


setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/bayesS_annot_chr")
setwd("/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/SBayesS")
plt_list = data.table()
for (chr in 1:22){
  res = fread(paste0("All_MVP_Trpchevska_De-Angelis_BBJ_filter_SbayesS_chr", chr, ".parRes"))
  res$chr = chr
  plt_list = rbind(plt_list, res[V1=="Pi" | V1=="S" | V1=="hsq", ])
}
colnames(plt_list)[1] = c("cate")
# fwrite(plt_list, "gctb_result.txt", sep="\t")

plt_list$y_max = plt_list$Mean + plt_list$SD
plt_list$y_min = plt_list$Mean - plt_list$SD
chr_len=fread("/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/chr_length.txt")
plt_dt = merge(plt_list, chr_len, by="chr")

h = sum(subset(plt_list, cate=="hsq")$Mean)
sm = mean(subset(plt_list, cate=="S")$Mean)
pim = mean(subset(plt_list, cate=="Pi")$Mean)

# pearson correaltion for hsq
h_dt = subset(plt_dt, cate=="hsq")
pearson_cor = cor(h_dt$Mean, h_dt$bp, method="pearson")
pearson_cor
[1] 0.9134548

pi = ggplot(subset(plt_list, cate=="Pi"), aes(x=chr, y=Mean)) +
  geom_bar(stat="identity", position="dodge", fill="#018B00") +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  scale_x_continuous(breaks = 1:22) +
  scale_y_continuous(expand=c(0, 0)) + 
  labs(x="CHR", y="Pi")+
  theme_classic()
ggsave("Pi_2.8m_annot_1M.pdf", pi, width=5, height=4)

hsq = ggplot(subset(plt_dt, cate=="hsq"), aes(x=bp, y=Mean)) +
  geom_smooth(aes(x=bp, y=Mean), method="lm", level=0.95)+
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  geom_point(aes(y = Mean), size = 0.5) +
  geom_text(aes(y = Mean, label = chr), hjust = -0.2, vjust = -0.5, size = 3) + 
  labs(x="Chromosome length (Mb)", y="Heritability")+
  theme_classic()
ggsave("hsq_2.8m_annot_1M.pdf", hsq, width=5, height=4)
 
s = ggplot(subset(plt_list, cate=="S"), aes(x=chr, y=Mean)) +
  geom_bar(stat="identity", position="dodge", fill="#B03060") +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.1) +
  scale_x_continuous(breaks = 1:22) +
  labs(x="CHR", y="S")+
  theme_classic()
ggsave("S_1M.pdf", s, width=5, height=4)