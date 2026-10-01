#! /bin/bash
RDIR="/public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript"
outdir="/public/share/wchirdzhq2022/Wulab_share/GWAS-summary/MVP/hearing-loss"

for file in *.GIA.dbGaP.txt.gz; do
    id=$(echo "${file}" | cut -d '.' -f 4)
    ${RDIR} reformat_MVP.r ${file} ${outdir}/MVP_${id}_reformat
    qsubshcom "$cmd" 1 100G reformat 90:00:00 ""
done

