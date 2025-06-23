#===========================================#
#                                           #
# use opera to generate the omics plot file #
#                                           #
#===========================================#
opera="/public/home/wuyang/program/OPERA-main/opera_Linux"
bfile="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/g1000_eur/g1000_eur"
mQTL="/public/home/gaikai/Age/data/mQTL/LBC_BSGS_meta/bl_mqtl_chr"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
Genelist="/public/home/shilulu/script/plot_smr/glist_hg19_refseq.txt"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"


# 17 ENSG00000004975 DVL2
# 17 ENSG00000170291 ELP5
# 19 ENSG00000178093 TSSK6
# 19 ENSG00000186010 NDUFA13
# 19 ENSG00000250067 YJEFN3
# 1 ENSG00000230896 RP11-767N6.7
# 6 ENSG00000220614 RP11-480N24.4

cat plot_probe | while read -r id; do
  arr=($id)
  chrs=${arr[0]}
  probe=${arr[1]}
  gene=${arr[2]}
  echo "${mQTL}${chrs}" > ${gene}_mlist
  cmd="$opera --bfile $bfile \
    --gwas-summary $gwas \
    --beqtl-summary $eQTLGen \
    --besd-flist ${gene}_mlist \
    --plot \
    --probe ${probe} \
    --probe-wind 2000 \
    --gene-list $Genelist \
    --out ${gene}_plot > ${gene}.opera.log 2>&1"
  qsubshcom "$cmd" 10 20G opera 90:00:00 ""
done

source("E:/Shi/ALL_of_my_Job/WCH_script/SMRplot/plot_OmicsSMR_xQTL.r")
setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/SMR")
SMRData = ReadomicSMRData("ACADVL_plot.ENSG00000072778.txt")
pdf("ACADVL_opera.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.19e-06, msmr_thresh=5.37e-07, m2esmr_thresh=5.75e-04, m2esmr_heidi=0.01,
                window=200, anno_methyl=TRUE, highlight="rs507506", annoSig_only=TRUE, max_anno_probe=6,
                eprobeNEARBY="ENSG00000072778",mprobeNEARBY=c("cg03613822", "cg12805420"), 
                epi_plot=TRUE, funcAnnoFile="E:/Shi/ALL_of_my_Job/WCH_script/SMRplot/funcAnno.RData")
dev.off()
# epi_plot=TRUE if you want plot epigenome
# anno_methyl=TRUE if you want the mQTL can annote in plot
pdf("ACADVL_DVL2_ELP5.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.19e-06, msmr_thresh=5.37e-07, m2esmr_thresh=5.75e-04, m2esmr_heidi=0.01,
                window=200, anno_methyl=TRUE, annoSig_only=TRUE, max_anno_probe=6,
                eprobeNEARBY=c("ENSG00000072778","ENSG00000170291","ENSG00000004975"),mprobeNEARBY=c("cg12805420"), 
                )
dev.off()




#------------- plotMhnSMR ----//
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)

dt = fread("ARHL_MVP_AJHG_BBJ_reformatMETAL_eQTLGen.smr")[, .(Gene, ProbeChr, Probe_bp, p_SMR, b_SMR, p_HEIDI)]
plt = dt[,.(Outco_Gene, Outco_Chr, Expo_bp, p_SMR)]
colnames(plt) = c("Gene", "CHR", "BP", "P")

Gene = dt[p_SMR < threshold & p_HEIDI > 0.01, ]
Gene$color = Gene[, ifelse(b_SMR > 0, brewer.pal(8,"Dark2")[1:2][2], brewer.pal(8,"Dark2")[1:2][1])]
Gene$pch = Gene[, ifelse(b_SMR > 0, 24, 25)]


plotMhnSMR <- function(plt, highlight=NULL, threshold=0.05/nrow(plt), ) {
  
  chr_lengths <- plt %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP))
  data <- plt %>%
    arrange(CHR, BP) %>%
    mutate(chr_cumsum = cumsum(c(0, head(chr_lengths$chr_len, -1)))[CHR]) %>%
    mutate(BPcum = BP + chr_cumsum)
  
  axis_set <- data %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  if (!is.null(highlight)){
    hlt = merge(highlight, data, by=("Gene"))
  }
  
}


hlt_dt = merge(Gene, data, by.x=c("Outco_Gene", "Expo_bp"), by.y=c("Gene", "BP"))
hlt_dt$pch = as.factor(hlt_dt$pch)

threshold = -log10(0.05 / nrow(dt))
p = ggplot() +
  geom_point(data = data, aes(x = BPcum, y = -log10(P), color = factor(CHR)), alpha = 0.9) +
  geom_point(data=hlt_dt, aes(x=BPcum, y=-log10(P), shape=pch),fill=hlt_dt$color, color=hlt_dt$color, alpha=0.9, size=2.5) +
  geom_text_repel(
    data = hlt_dt,
    aes(x=BPcum, y=-log10(P), label = Outco_Gene),
    color=hlt_dt$color,
    size = 3,
    fontface="italic",
    max.overlaps = Inf,
    max.iter = 10000, 
    force = 2,
    force_pull = 0.01,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50") + 
  geom_hline(yintercept=-log10(threshold), color="grey68", linetype="dashed") +
  scale_color_manual(values=rep(c("#bfbfbf", "#b3dfe7"), 22)) + 
  scale_shape_manual(values = c(24, 25)) + 
  scale_x_continuous(label=axis_set$CHR, breaks=axis_set$center, expand=expansion(mult = c(0, 0.01))) +
  scale_y_continuous(expand=expansion(mult = c(0, 0.01))) + 
  labs(x="Chromosome", y=expression(-log[10] (italic(P_SMR)))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=15, margin=margin(t=10), color="black"),
        axis.title.y = element_text(size=15, margin=margin(t=10), color="black"),
        axis.line = element_line(color="black"), 
        axis.ticks = element_line(color="black"))
p
ggsave("mQTL_mnh.png", p, height = 4, width = 10, dpi=600)
ggsave("mQTL_mnh.pdf", p, height = 4, width = 10)