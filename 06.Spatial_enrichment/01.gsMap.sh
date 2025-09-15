# spatial LDSC
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
gwas="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/07.spEnrich/All_MVP_Trpchevska_De-Angelis_BBJ_filter.sumstats.gz"
hm3_w="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/weights_hm3_no_hla/weights."
gsmap run_spatial_ldsc \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_sup' \
    --sumstats_file ${gwas} \
    --w_file ${hm3_w} \
    --num_processes 30

# cauchy combination
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
gsmap run_cauchy_combination \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_sup' \
    --annotation 'annotation'

# report result
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
gwas="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/07.spEnrich/All_MVP_Trpchevska_De-Angelis_BBJ_filter.sumstats.gz"
gsmap run_report \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL_sup' \
    --annotation 'annotation' \
    --sumstats_file ${gwas} \
    --top_corr_genes 50
