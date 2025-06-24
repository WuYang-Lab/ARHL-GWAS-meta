#*******************************************************************************#
#**  GCTB
# We introduced a method (called BayesS) to infer the action of natural selection 
# on the genetic variants underlying a complex trait. 
# stimating the relationship between the variance of SNP effects and MAF 
# Zeng et al. Nat Gent
#*******************************************************************************#
cat ukbEURu_hm3_sparse_mldm_list.txt | while read -r id; do
	dir="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEURu_hm3_shrunk_sparse"
	file=$(basename "$id")
	echo "${dir}/${file}" >> ukbEURu_hm3_sparse_mldm_list_my.txt
done


#************************************#
#sbayes estimate genetic archtecture #
#************************************#

# # # 1 step: imputed the summary data
#! /bin/bash
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL.txt"
out_block="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/impu_block/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_impu"

# step1: impute summary data --block
cmd1="${gctb} --ldm-eigen ${ldref} --gwas-summary ${gwas} --impute-summary --block {TASK_ID} --out ${out_block} --thread 10 > ${out_block}{TASK_ID}.log 2>&1"
qsubshcom "$cmd1" 10 100G gctb_impu 90:00:00 "-array=1-591"

out_impu="ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb.imputed"
# step2: combined the all block file
cd ./impu_block
gctb --gwas-summary ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_impu --merge-block-gwas-summary --out ${out_impu} > ../${out_impu}.log 2>&1

mv ${out_impu}.ma ../${out_impu}.ma

# # # 2 step: run sBayeS
# # # 2.8 M snp + baseline annot
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukb_50k_bigset_2.8M"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/bayesS_annot_chr"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
gwas_impu="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb.imputed.ma"
outprx=$(basename -- ${gwas_impu} ".imputed.ma")

cmd="gctb --sbayes S --ldm ${ldref}/ukb50k_shrunk_chr{TASK_ID}_mafpt01.ldm.sparse --annot ${annofile} --chr {TASK_ID} \
--gwas-summary ${gwas_impu} --no-mcmc-bin --thread 10 \
--out ${outdir}/${outprx}_SbayesS_chr{TASK_ID} > ${outdir}/${outprx}_SbayesS_chr{TASK_ID}.log 2>&1"
qsubshcom "$cmd" 10 100G BayesS 90:00:00 "-array=1-22"