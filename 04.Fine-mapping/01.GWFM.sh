#----------------------------------------------------#
# # # 1 step: imputed the summary data
#! /bin/bash
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
gwas="/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt"
out_block="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/impu_block/All_MVP_Trpchevska_De-Angelis_BBJ_impu"

# step1: impute summary data --block split block running is a efficiency way, totally 591 block 
for block in {1..591}; do
${gctb} --ldm-eigen ${ldref} \
--gwas-summary ${gwas} \
--impute-summary \
--block ${block} \
--thread 10 \
--out ${out_block} > ${out_block}${block}.log 2>&1
done

out_impu="All_MVP_Trpchevska_De-Angelis_BBJ.imputed"
# step2: combined the all block file to one 
cd ./impu_block
$gctb --gwas-summary All_MVP_Trpchevska_De-Angelis_BBJ_impu --merge-block-gwas-summary --out ${out_impu} > ../${out_impu}.log 2>&1


# Run genome-wide fine-mapping analysis
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
outdir="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/gwfm"
gwas_impu="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/All_MVP_Trpchevska_De-Angelis_BBJ.imputed.ma"

outprx=$(basename -- ${gwas_impu} ".imputed.ma")
${gctb} --sbayes RC \
--ldm-eigen ${ldref} \
--annot ${annofile} \
--gwas-summary ${gwas_impu} \
--n-dist-auto \
--write-mcmc-bin \
--thread 10 \
--out ${outdir}/${outprx}_SbayesRC > ${outdir}/${outprx}_SbayesRC.log 2>&1

# # # Calculate credible sets
# step2: use ldfile calculate credible sets
gctb="/public/home/shilulu/software/gctb"
ldfile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukbEUR_Imputed/LD_rsq05.ld.txt"
mcmc_prx="All_MVP_Trpchevska_De-Angelis_BBJ_SbayesRC"

$gctb --cs \
--ld-file ${ldfile} \
--pip 0.9 \
--mcmc-samples ${mcmc_prx} \
--out ${mcmc_prx}_finemapping > ${mcmc_prx}_finemapping.log 2>&1