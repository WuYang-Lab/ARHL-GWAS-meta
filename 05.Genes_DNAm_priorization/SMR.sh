#**********************#
# SMR & HEIDI analysis #
#**********************#
# whole blood eQTL and mqtl#
smr="/public/home/shilulu/software/smr-1.3.1/smr"
bfile="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
mQTL="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/mQTL/LBC_BSGS_meta_all"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
gwas="/public/home/shilulu/Wulab/sll/ARHL/NC_sup_test/01.meta/All_MVP_Trpchevska_De-Angelis_BBJ_filter.txt"
wkdir="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/05.SMR"
prefix=$(basename "$gwas" ".txt")

qsubshcom "$smr --bfile $bfile \
--beqtl-summary $eQTLGen \
--gwas-summary $gwas \
--diff-freq-prop 0.5 \
--maf 0.01 \
--heidi-min-m 3 \
--heidi-max-m 10 \
--cis-wind 2000 \
--thread-num 10 \
--out $wkdir/eqtlgen/$prefix \
> $wkdir/eqtlgen/${prefix}_eQTLGen.log 2>&1" 10 100G smr_eqtlgen 1:00:00 ""



#=====================================================================================#
#                                                                                     #
# two molecular traits test of whole blood methQTL and eQTL, only for significant QTL #
#                                                                                     #
#=====================================================================================#
# --beqtl-summary the first one reads mQTL summary data as the exposure. 
# The second one reads eQTL summary data from as the outcome.

# run smr of mqtl to eqtl 
#! /bin/bash 
smr="/public/home/shilulu/software/smr-1.3.1/smr"
bfile="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
mQTL="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/mQTL/LBC_BSGS_meta_all"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
mprobe="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/05.SMR/mqtl/mQTL_sign_probeID.txt"
eprobe="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/05.SMR/eqtlgen/eQTL_sign_probeID.txt"
outdir="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/05.SMR/m2eqtl"

qsubshcom "$smr --bfile ${bfile} \
--beqtl-summary ${mQTL} \
--beqtl-summary ${eQTLGen} \
--extract-exposure-probe $mprobe \
--extract-outcome-probe $eprobe \
--diff-freq-prop 0.5 \
--maf 0.01 \
--heidi-min-m 3 \
--heidi-max-m 10 \
--cis-wind 2000 \
--thread-num 10 \
--out ${outdir}/m2eqtl \
> ${outdir}/m2eqtl.log 2>&1" 1 100G m2eqtl 1:00:00 ""