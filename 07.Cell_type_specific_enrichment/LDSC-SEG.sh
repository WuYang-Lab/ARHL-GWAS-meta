conda activate ldsc
#*************************#
##       LDSC-SEG        ##
#*************************#
#! /bin/bash
cd /public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/03.ldsc
sums="~/software/ldsc/munge_sumstats.py"
SNPlist="~/Wulab/LDSC/eur_w_ld_chr/w_hm3.snplist"
gwas="~/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.gz"
prefix=$(basename $gwas ".gz")

$sums --sumstats $gwas \
--merge-alleles $SNPlist \
--chunksize 500000 \
--a1 A1 \
--a2 A2 \
--out $outdir/$prefix


#***************************************************************#
## LDSC-SEG for Mouse cochleae in Philippe Jean et al. PNAS
#***************************************************************#
conda activate ldsc

Eshel="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Eshel/LDSC/config/scRNA_gene_expr.ldcts"
Milon="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Milon/LDSC/config/scRNA_gene_expr.ldcts"
Ranum="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Ranum/LDSC/config/scRNA_gene_expr.ldcts"
Jean="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Jean/LDSC/config/scRNA_gene_expr.ldcts"
Sun1_2="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Sun/1-2M/LDSC/config/scRNA_gene_expr.ldcts"
Sun12_15="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Sun/12-15M/LDSC/config/scRNA_gene_expr.ldcts"
Sun5="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Sun/5M/LDSC/config/scRNA_gene_expr.ldcts"
Xu="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Xu/LDSC/config/scRNA_gene_expr.ldcts"
Iyer="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/08.scEnrich/Iyer/LDSC/config/scRNA_gene_expr.ldcts"

REF_LD_CTS=$Iyer
REF_LD_CHR="~/Wulab/LDSC/1000G_EUR_Phase3_baseline/baseline."
W_LD_CHR="~/Wulab/LDSC/weights_hm3_no_hla/weights."
gwas="~/Wulab/sll/ARHL/NC_sup_test/03.ldsc/All_MVP_Trpchevska_De-Angelis_BBJ_filter.sumstats.gz"
$ldsc --h2-cts $gwas \
        --ref-ld-chr $REF_LD_CHR \
        --out All_MVP_Trpchevska_De-Angelis_BBJ_filter \
        --ref-ld-chr-cts $REF_LD_CTS \
        --w-ld-chr $W_LD_CHR

