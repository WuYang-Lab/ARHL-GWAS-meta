#************************************#
#sbayes estimate genetic archtecture #
#************************************#
gctb="/public/home/shilulu/software/gctb"
ldref="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/ukb_50k_bigset_2.8M"
annofile="/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt"
gwas="/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt"
outdir="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/04.gctb/SBayesS"
outprx=$(basename -- ${gwas} ".txt")

for chr in {1..22}; do
$gctb --sbayes S \
--ldm ${ldref}/ukb50k_shrunk_chr${chr}_mafpt01.ldm.sparse \
--annot ${annofile} \
--chr ${chr} \
--gwas-summary ${gwas} \
--no-mcmc-bin \
--thread 10 \
--out ${outdir}/${outprx}_SbayesS_chr${chr} > ${outdir}/${outprx}_SbayesS_chr${chr}.log 2>&1
done