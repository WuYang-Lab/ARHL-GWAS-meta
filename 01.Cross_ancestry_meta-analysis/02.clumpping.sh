###########################
## PLINK1.9  LD clumping ##
# hearing loss
###########################
# using eur reference genotype
eurdir="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur"
plink \
    --clump All_MVP_Trpchevska_De-Angelis_BBJ_filter.gz \
    --bfile ${eurdir}/g1000_eur \
    --clump-p1 5e-8 \
    --clump-p2 5e-8 \
    --clump-r2 0.05 \
    --clump-kb 1000 \
    --clump-snp-field SNP \
    --clump-field p \
    --out ./clump/All_MVP_Trpchevska_De-Angelis_BBJ_filter.clump
