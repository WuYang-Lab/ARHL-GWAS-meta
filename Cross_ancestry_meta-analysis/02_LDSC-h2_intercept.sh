# LDSC intercept
#*************************#
## estimate heritability ##
#*************************#
conda activate ldsc

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

# step 2: calculate h2 
LDSC="/public/home/shilulu/software/ldsc/ldsc.py"
REF_LD_CHR="/public/home/lijunpeng/LDSC_LD_data/eur_w_ld_chr/"

cmd2="${LDSC} --h2 ${outdir}/ARHL_MVP_AJHG_BBJ_ldsc.sumstats.gz \
        --ref-ld-chr ${REF_LD_CHR} \
        --w-ld-chr ${REF_LD_CHR} \
        --out ${outdir}/ARHL_MVP_AJHG_BBJ_h2_observe"
qsubshcom "$cmd2" 1 100G ldsc_h2 90:00:00 "-wait=${sub_id1}" 

#***********************************#
## convert to liability for binary ##
# only need to run this rather above
#***********************************#
# data     N     case   control
# total	1474404	438594	1035810

LDSC="/public/home/shilulu/software/ldsc/ldsc.py"
REF_LD_CHR="/public/home/lijunpeng/LDSC_LD_data/eur_w_ld_chr/"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc"

cmd2="${LDSC} \
--h2 ${outdir}/ARHL_MVP_AJHG_BBJ_ldsc.sumstats.gz \
--ref-ld-chr ${REF_LD_CHR} \
--w-ld-chr ${REF_LD_CHR} \
--out ${outdir}/ARHL_MVP_AJHG_BBJ_h2_0.4 \
--samp-prev 0.297 \
--pop-prev 0.05"
qsubshcom "$cmd2" 1 100G ldsc_h2 90:00:00 "-wait=${sub_id1}" 