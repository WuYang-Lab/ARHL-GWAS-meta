opera="/public/home/wuyang/program/OPERA-main/opera_Linux"
bfile="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
mQTL="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/mQTL/LBC_BSGS_meta_all"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
Genelist="/public/home/shilulu/script/plot_smr/glist_hg19_refseq.txt"
gwas="/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt"

echo "$mQTL" > mlist
# ENSG00000072778 ACADVL
# ENSG00000004975 DVL2
# ENSG00000170291 ELP5
declare -A probe=(
  ["ACADVL"]="ENSG00000072778"
  ["DVL2"]="ENSG00000004975"
  ["ELP5"]="ENSG00000170291"
)

for gene in "${!probe[@]}"; do
  en=${probe[$gene]}
  echo "create plot probe ${en} for ${gene}"
  $opera --bfile $bfile \
  --gwas-summary $gwas \
  --beqtl-summary $eQTLGen \
  --besd-flist mlist \
  --plot \
  --probe $en \
  --probe-wind 500 \
  --gene-list $Genelist \
  --out ${gene}_plot > ${gene}.opera.log 2>&1
done


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