source("F:/Shi/ALL_of_my_Job/WCH_script/plot/SMRplot/plot_OmicsSMR_xQTL.r")
setwd("F:/Shi/ALL_of_my_Job/24-28 Ph.D WCHSCU/2_project_hearing loss/NC_revision/SMR")
SMRData = ReadomicSMRData("ACADVL_plot.ENSG00000072778.txt")
pdf("ACADVL.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.20e-06, msmr_thresh=5.39e-07, m2esmr_thresh=5.38e-04, m2esmr_heidi=0.01,
                 window=200, anno_methyl=TRUE, annoSig_only=TRUE, max_anno_probe=8,
                 eprobeNEARBY="ENSG00000072778",mprobeNEARBY=c("cg00072720", "cg12805420"),
                 epi_plot=TRUE, funcAnnoFile="F:/Shi/ALL_of_my_Job/WCH_script/plot/SMRplot/funcAnno.RData")
dev.off()

# epi_plot=TRUE if you want plot epigenome
# anno_methyl=TRUE if you want the mQTL can annote in plot
pdf("ACADVL_DVL2_ELP5.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.20e-06, msmr_thresh=5.39e-07, m2esmr_thresh=5.38e-04, m2esmr_heidi=0.01,
                window=200, anno_methyl=TRUE, annoSig_only=TRUE, max_anno_probe=8,
                eprobeNEARBY=c("ENSG00000072778","ENSG00000170291","ENSG00000004975"),mprobeNEARBY=c("cg12805420"))
dev.off()