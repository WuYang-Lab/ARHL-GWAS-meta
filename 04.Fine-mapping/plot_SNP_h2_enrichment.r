# plt 
library(data.table)
setwd("F:/Github/PHD_job/2_project_hearing loss/NC_revision/gctb/gwfm")
res = fread("All_MVP_Trpchevska_De-Angelis_BBJ_SbayesRC.parSetRes",head=F,check.names=F)
plotRes = subset(res, V1 == "AnnoPerSnpHsq_Enrichment")
plotRes1 = subset(res, V1 == "AnnoJointProb_pi1")
plotRes2 = subset(res, V1 == "AnnoJointProb_pi2")
plotRes3 = subset(res, V1 == "AnnoJointProb_pi3")
plotRes4 = subset(res, V1 == "AnnoJointProb_pi4")
plotRes5 = subset(res, V1 == "AnnoJointProb_pi5")

bottom5 = sort(plotRes$V3)[5]; top20 = sort(plotRes$V3, decreasing=T)[20] 
idx = which(plotRes$V3 >= top20 | plotRes$V3 <= bottom5)

#   combine plot data
plotData = data.frame(rbind(plotRes[idx,c(2,3,4)], plotRes1[idx,c(2,3,4)], plotRes2[idx,c(2,3,4)], 
                            plotRes3[idx,c(2,3,4)],plotRes4[idx,c(2,3,4)], plotRes5[idx,c(2,3,4)]))

plotData$compLabel = c(rep("Hsq enrich", length(idx)), rep("Zero effect", length(idx)), rep("Small effect", length(idx)),
                       rep("Medium effect", length(idx)), rep("Large effect", length(idx)), rep("Very large effect", length(idx)))

colnames(plotData)[1:3] = c("cates","enrich","enrich_se")
plotData$y_min = ifelse(plotData$enrich > 0, plotData$enrich, plotData$enrich - plotData$enrich_se)
plotData$y_max = ifelse(plotData$enrich > 0, plotData$enrich + plotData$enrich_se, plotData$enrich)

catesNew = c("Conding","Conserved_LindbladToh","Enhancer_Andersson","H3K4me3",
             "H3K9ac","Intron","Repressed","SuperEnhancer","Transcription","TSS",
             "UTR_3","UTR_5","GERP","MAFbin3","MAFbin8","Synonymous","Non_synonymous","Conversed_Vertebrate","Conserved_Mammal","Conserved_Primate",
             "BivFlank","Promoter","Promoter_ancient","Promoter_ancient.flanking","Promoter_ExAC")

plotData$cates=catesNew;
plotData$cates = factor(plotData$cates, levels = plotData$cates[order(plotData$enrich[1:length(idx)], decreasing = FALSE)])
plotData$compLabel = factor(plotData$compLabel, levels=unique(plotData$compLabel))

require(ggplot2); require(RColorBrewer); require(dplyr)

rder_label = plotData$cates[order(plotData$enrich[1:length(idx)], decreasing = FALSE)]
plotData$cates = factor(plotData$cates, levels = rder_label)

levels(plotData$compLabel)[1] = "Per-SNP heritability enrichment"
my_color1 = c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(5, "Set2"))
my_color = rev(my_color1)



subplot = subset(plotData, compLabel!="Per-SNP heritability enrichment")
target_compLabel = "Zero effect"
baseline = 0.98
subplot = subplot %>%
  mutate(
    enrich_adj = ifelse(compLabel == target_compLabel, enrich - baseline, enrich),
    y_min_adj = ifelse(compLabel == target_compLabel, y_min - baseline, y_min),
    y_max_adj = ifelse(compLabel == target_compLabel, y_max - baseline, y_max),
    y_min_adj = pmax(y_min_adj, 0),
    y_max_adj = pmin(y_max_adj, ifelse(compLabel == target_compLabel, 1 - baseline, Inf)))



common_theme = theme(axis.text.x = element_text(size = 12, angle=60, hjust= .5, vjust= 0.5),
                     axis.text.y = element_text(size = 12, angle=0),
                     plot.title = element_text(hjust = 0.5),
                     axis.line = element_line(colour = "black"),
                     strip.background = element_rect(colour = NA),
                     strip.text = element_text(size = 9))

# 获取所有 compLabel 值（假设有 5 个）
comp_labels = unique(subplot$compLabel)
if (length(comp_labels) != 5) warning("Expected 5 compLabels, found ", length(comp_labels))

# 创建子图列表
plot_list = list()
for (label in comp_labels) {
  sub_data = subplot %>% filter(compLabel == label)
  if (label == target_compLabel) {
    # Zero effect 子图
    p = ggplot(sub_data, aes(x = cates, y = enrich_adj, fill = cates)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = y_min_adj, ymax = y_max_adj, color = cates)) +
      guides(fill = "none", color = "none") +
      labs(x = "Annotation", y = "", title = label) +
      coord_flip() +
      theme_bw()+
      common_theme +
      scale_y_continuous(
        limits = c(0, 0.02),
        breaks = seq(0, 0.02, by = 0.005),
        labels = sprintf("%.3f", seq(0.98, 1, by = 0.005)))
  } else {
    # 其他子图
    p = ggplot(sub_data, aes(x = cates, y = enrich, fill = cates)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = y_min, ymax = y_max, color = cates)) +
      guides(fill = "none", color = "none") +
      labs(x="", y="", title = label) +
      coord_flip() +
      theme_bw()+
      common_theme+
      theme(axis.text.y = element_blank())
  }
  plot_list[[label]] = p
}

# 使用 wrap_plots 组合 5 个子图，固定宽度
combined_plot = wrap_plots(plot_list, nrow = 1, widths = rep(1, length(plot_list)))

# 输出图形
print(combined_plot)
ggsave("SBayesRC_hsq_enrich.pdf", combined_plot, height=8, width=15)



# per snp enrichment
per_snp = subset(plotData, compLabel=="Per-SNP heritability enrichment")
library(ggplot2)
p1 = ggplot(per_snp, aes(x=cates, y=enrich, fill=cates)) +
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=y_min, ymax=y_max, color=cates)) + 
  scale_y_continuous(expand=c(0,0)) +
  guides(fill="none", color="none") +  
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="Annotation", y="Per SNP heritability enrichment") +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_text(size=10, angle=0),
        axis.text.y = element_text(size=10, angle=0))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        strip.background = element_rect(colour=NA))
p1
ggsave("per_SNP_enrichment.pdf", p1, height=7, width=7)
