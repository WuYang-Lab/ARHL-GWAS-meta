#***********************************#
## LDSC-SEG ##
#***********************************#
conda activate ldsc

# # # # # # # # # # # # # # # # # # # # # # # # 

# genetrate own annotate file and run ldsc

# # # # # # # # # # # # # # # # # # # # # # # # 
## Step 1: Creating an annot file
# # # 1) generate matrix in python
conda activate gsMap
python -c ""
import pandas as pd 
import scanpy as sc 

SC = sc.read_h5ad("ScRNA-seq-P8_P12_P20_mouse_cochlea_Jun27.h5ad")
expr_df = SC.to_df()
expr_df['cell_type'] = SC.obs['cell_type'].values

# per gene avrg expression in per cell type, rm expr 0 in all cell type
mean_expr_in_each_cell_type = expr_df.groupby('cell_type').mean().T
mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]
mean_expr_in_each_cell_type.to_csv('Sc_P8_12_20.txt', sep="\t")

# for different postnatal ages
expr_df['age'] = SC.obs['age'].values
# per gene avrg expression in per cell type, rm expr 0 in all cell type
age_8 = expr_df[expr_df['age'] == "P8"]
expr_age_8 = age_8.drop(columns=['age'])
age_12 = expr_df[expr_df['age'] == "P12"]
expr_age_12 = age_12.drop(columns=['age'])
age_20 = expr_df[expr_df['age'] == "P20"]
expr_age_20 = age_20.drop(columns=['age'])

mean_expr_in_each_cell_type = expr_age_20.groupby('cell_type').mean().T
mean_expr_in_each_cell_type = mean_expr_in_each_cell_type.loc[(mean_expr_in_each_cell_type != 0).any(axis=1)]
mean_expr_in_each_cell_type.to_csv('P20/Sc_P20.txt', sep="\t")
#================================================================

# # # 2) generate homologous gene in R
# if (!require("BiocManager", quietly = TRUE)) 
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(dplyr)
library(data.table)
library(biomaRt)
listEnsembl()
# connect Ensembl
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
## 构建转化函数
homologs <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name",
                                "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"),
                  filters = "with_mmusculus_homolog", 
                  values = TRUE, 
                  mart = human)
# head(homologs)
homologs_1to1 <- subset(homologs, mmusculus_homolog_orthology_type == "ortholog_one2one")
# rm MT X Y chromosome genes
homGene_1to1_chr1_22 <- homologs_1to1 %>% 
  filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
  rename(human_ensembl_gene = ensembl_gene_id)

fwrite(homGene_1to1_chr1_22, file="HomGene1to1_mouse2human.txt", sep="\t", row.names=FALSE)
#================================================================

# # # 3) extract 1:1 homologs scRNA genes from mouse cochlea
# 同时生成control.GeneSet文件，也就是所有的基因
library(dplyr)
library(data.table)
setwd("/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/")
scGene = fread("Sc_P8_12_20.txt")
homologs_1to1 = fread("HomGene1to1_mouse2human.txt")[,2:3,with=FALSE]
colnames(homologs_1to1) = c("human_gene", "gene")

merge_df = merge(scGene, homologs_1to1, by="gene")
merge_df = subset(merge_df, select = -c(gene))
colnames(merge_df) = gsub("[[:space:]/-]", "_", colnames(merge_df))
colnames(merge_df) = gsub("'", "_", colnames(merge_df))
merge_df = merge_df %>% select(human_gene, everything())

control = merge_df[, 1]
fwrite(control, file="Sc_P8_12_20_control.GeneSet", sep="\t", col.names=FALSE)

lapply(2:ncol(merge_df), function(i){
  temp_dt = data.frame(merge_df[[1]], merge_df[[i]])
  # colnames(df1) = names(merge_df)[c(37, i)]
  sorted_dt = temp_dt[order(temp_dt[[2]], decreasing = TRUE), ]
  top_10 = sorted_dt[1:floor(0.1*nrow(sorted_dt)), ]
  top_10_gene = top_10[, 1, drop=FALSE]

  cell_type = colnames(merge_df)[i]
  fwrite(top_10_gene, file=paste0("Sc_P8_12_20_", cell_type, ".GeneSet"), sep="\t", col.names=FALSE)
})
#================================================================

# # # 4) use make_annot.py to generate annot.gz file 
# using 100K widow as Hilary K. Finucane et al. Nat Genet

# generate ENSG.coord.txt
library(biomaRt)
library(data.table)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
human_gene <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name"),
                    filters = "", values = TRUE, mart = human)
colnames(human_gene) = c("CHR", "ENSG_ID", "GENE")

glist_hg19 <- fread("/public/home/shilulu/script/plot_smr/glist-hg19")
colnames(glist_hg19) = c("CHR", "START", "END", "GENE")
merge_df <- merge(glist_hg19, human_gene, by="GENE")
glist_hg19_ENSG <- merge_df[, c("ENSG_ID", "CHR.y", "START", "END")]
colnames(glist_hg19_ENSG) = c("GENE", "CHR", "START", "END")

fwrite(glist_hg19_ENSG, file="ENSG_coord.txt", sep="\t")

# run make_annot.py 为control 和 每个cell type生成
conda activate ldsc
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  # for chr in {1..22}; do 
  cmd="python ${ldsc}/make_annot.py \
      --gene-set-file ${annot}.GeneSet \
      --gene-coord-file ${ldsc}/ENSG_coord.txt \
      --bimfile ${bfile}.{TASK_ID}.bim \
      --windowsize 100000 \
      --annot-file ${annot}.{TASK_ID}.annot.gz"
  # done 
  qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22"
done

## Step 2: Computing LD scores with an annot file
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/share/wchirdzhq2022/Wulab_share/LDSC"

awk '{if ($1!="SNP") {print $1} }' ${ldsc_dir}/w_hm3.snplist > ${ldsc_dir}/listHM3.txt

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  cmd="python ${ldsc}/ldsc.py --l2 \
    --bfile ${bfile}.{TASK_ID} \
    --print-snps ${ldsc_dir}/listHM3.txt \
    --ld-wind-cm 1 \
    --annot ${annot}.{TASK_ID}.annot.gz \
    --thin-annot \
    --out ${annot}.{TASK_ID}"
  qsubshcom "$cmd" 1 10G ldsc_l2 1:00:00 "-array=1-22"
done

# generate .ldcts file for my sc data
dir=`pwd`
ls *GeneSet | while read id; do
  tissue=$(basename -- ${id} | sed 's/^Sc_P20_//;s/\.GeneSet$//')
  if [[ "$tissue" == "control" ]]; then
    echo "control not to be the tissue or cell type"
  else
    cell=${dir}/$(basename -- ${id} "GeneSet")
    control=${dir}/"Sc_P20_control."
    echo "${tissue} ${cell},${control}" >> Mouse_cochleae_P20_gene_expr.ldcts
  fi
done
# rm control file in cell 

#***************************************************************#
## LDSC-SEG for Mouse cochleae in Philippe Jean et al. PNAS
#***************************************************************#
#! /bin/bash
SUMS="/public/home/shilulu/software/ldsc/munge_sumstats.py"
SNPlist="/public/home/lijunpeng/LDSC_LD_data/eur_w_ld_chr/w_hm3.snplist"
meta_file="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL.gz"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc"

# step 1: covert to .sumstats.gz
cmd1="${SUMS} --sumstats ${meta_file} \
--merge-alleles ${SNPlist} \
--chunksize 500000 \
--a1 A1 \
--a2 A2 \
--out ${outdir}/ARHL_MVP_AJHG_BBJ_ldsc"
sub_id1=`qsubshcom "$cmd1" 1 100G sums_ldsc 90:00:00 ""`


#! /bin/bash                     
LDSC="/public/home/shilulu/software/ldsc/ldsc.py"
REF_LD_CHR="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_baseline/baseline."
REF_LD_CTS="/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/Mouse_cochleae_gene_expr.ldcts"
W_LD_CHR="/public/share/wchirdzhq2022/Wulab_share/LDSC/weights_hm3_no_hla/weights."
out_sums="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc"
out_ldscseg="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc/Mouse_cochleae"

cmd="${LDSC} --h2-cts ${out_sums}/ARHL_MVP_AJHG_BBJ_ldsc.sumstats.gz \
             --ref-ld-chr ${REF_LD_CHR} \
             --out ${out_ldscseg}/ARHL_MVP_AJHG_BBJ \
             --ref-ld-chr-cts ${REF_LD_CTS} \
             --w-ld-chr ${W_LD_CHR}"
qsubshcom "$cmd" 1 50G ldsc_cts 10:00:00 ""


#===============================================================#
####     VLOOKUP in R and calculate FDR and plot LDSC cts    ####
#===============================================================#
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript

data <- read.table('ARHL_MVP_AJHG_BBJ.cell_type_results.txt', header=TRUE, sep='\t')

data$FDR <- p.adjust(data$Coefficient_P_value, method='BH')
# save the FDR file
write.table(data, file='ARHL_MVP_AJHG_BBJ.cell_type_results.FDR.txt', sep='\t', row.name=FALSE, quote=FALSE)

library(ggplot2)
library(data.table)

data <- fread("ARHL_MVP_AJHG_BBJ.cell_type_results.FDR_plot.txt")
if (any(data$FDR <= 0.05)) {
  # 计算FDR <= 0.05的值中，对应的max p值，作为 fig 中的阈值线
  p_value_threshold <- max(data$Coefficient_P_value[data$FDR <= 0.05], na.rm=TRUE)
  
  log10_p_value_threshold <- -log10(p_value_threshold)
} else {
  log10_p_value_threshold <- NaN
}
# color for dif tissue
custom_colors <- c("Hair cells" = "#DD7694", 
                   "Supporting cells" = "#BCD4E7", 
                   "Surrounding structures" = "#056E83", 
                   "Lateral wall" = "#E9D9BF", 
                   "Circulating cells" = "#D4920A", 
                   "Glial cells" = "#5AA4AE", 
                   "Neurons" = "#65472F")   
# set Name as factor
data$Name <- factor(data$Name, levels=unique(data$Name))
data$Region <- factor(data$Region, levels=unique(data$Region))

p=ggplot(data, aes(x=Name, y=-log10(Coefficient_P_value), fill=Region)) +
  geom_bar(stat="identity", position="dodge", width=0.8)+
  geom_hline(yintercept=log10_p_value_threshold, linetype="dashed", color="red")+
  scale_fill_manual(values=custom_colors) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3))+
  labs(x=NULL, y="-log10(P)", title="Cells from Mouse cochlear")+
  theme_bw()+
  theme(axis.text.x=element_text(size=9, angle=60, hjust=1, vjust=1),
        axis.text.y=element_text(size=9),
        legend.position="left",
        legend.title=element_blank(),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))

ggsave("Mouse_cochlear_ldsc-seg.png", p, width=8, height=4, dpi=500)
