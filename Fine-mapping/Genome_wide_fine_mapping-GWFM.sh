#******************************#
#                              #
#   sBayes RC fine mapping     #
#                              #
#******************************#

# # # 1 step: imputed the summary data
#! /bin/bash
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL.txt"
out_block="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/impu_block/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_impu"

# step1: impute summary data --block 分blcok跑，总共是591个block，效率更高
cmd1="${gctb} --ldm-eigen ${ldref} --gwas-summary ${gwas} --impute-summary --block {TASK_ID} --out ${out_block} --thread 10 > ${out_block}{TASK_ID}.log 2>&1"
qsubshcom "$cmd1" 10 100G gctb_impu 90:00:00 "-array=1-591"

out_impu="ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb.imputed"
# step2: combined the all block file
cd ./impu_block
gctb --gwas-summary ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_impu --merge-block-gwas-summary --out ${out_impu} > ../${out_impu}.log 2>&1

mv ${out_impu}.ma ../${out_impu}.ma



# # # 2 step: run genome-wide fine-mapping analysis
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping"
gwas_impu="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb.imputed.ma"


outprx=$(basename -- ${gwas_impu} ".imputed.ma")
cmd="${gctb} --sbayes RC --ldm-eigen ${ldref} --annot ${annofile} --gwas-summary ${gwas_impu} --n-dist-auto --write-mcmc-bin --thread 10 \
--out ${outdir}/${outprx}_SbayesRC > ${outdir}/${outprx}_SbayesRC.log 2>&1"
qsubshcom "$cmd" 10 100G sbayesRC 90:00:00 ""

ldfile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed/LD_rsq05.ld.txt"
mcmc_prx=${outdir}/${outprx}_SbayesRC

cmd="gctb --cs --ld-file ${ldfile} --pip 0.9 --mcmc-samples ${mcmc_prx} --out ${mcmc_prx}_finemapping > ${mcmc_prx}_finemapping.log 2>&1"
qsubshcom "$cmd" 10 100G sbayesRC 90:00:00 ""

#=============================================================#
# # # generate SNP annotation of  Number of 1-SNP credible sets
#=============================================================#
# read local credible set
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/")
lcs = "ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC_finemapping.lcs"
library(data.table)
library(dplyr)

lcs = fread('ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC_finemapping.lcs')[, .(size, SNP)][, SNP := gsub("[\t\r\n]", "", SNP)]; 
lcs_snp = lcs[size==1, SNP]
baseline2 = fread('/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt')
lcs_snp_annot = baseline2[SNP %in% lcs_snp]
fwrite(lcs_snp_annot, 'finemapping_snp.txt', col.names=TRUE, sep='\t')