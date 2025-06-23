#===========================================================================#
#                               *  gsMap  *                                 #
#                             *             *                               #
#                           *                 *                             #
#      Genetically informed spatial mapping of cells for complex traits     #
#                                                                           #
#    Liyang Song, Wenhao Chen, Junren Hou, Minmin Guo, Jian Yang (2024)     #
#  Spatially resolved mapping of cells associated with human complex traits #
#===========================================================================#
conda activate gsMap

# 01. data format
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap"
outprx=$(basename -- "$gwas")
gsmap format_sumstats --sumstats ${gwas} --snp SNP --a1 A1 --a2 A2 --frq freq --beta beta --se SE --p p --n N \
--out ${outdir}/${outprx} > ${outdir}/${outprx}.log 2>&1

# 1. find latent representations
# The --workdir parameter specifies the working directory for gsMap, where all outputs will be saved.
ST_name="E16.5_E2S11.MOSTA"
ST_sample="/public/share/wchirdzhq2022/Wulab_share/gsMap/gsMap_example_data/ST/E16.5_E2S11.MOSTA.h5ad"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
gsmap run_find_latent_representations \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --input_hdf5_path ${ST_sample} \
    --annotation 'annotation' \
    --data_layer 'count'

# 2. generate gene specificity scores
homo="/public/share/wchirdzhq2022/Wulab_share/gsMap/homologs/mouse_human_homologs.txt"
ST_name="E16.5_E2S11.MOSTA"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
gsmap run_latent_to_gene \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --annotation 'annotation' \
    --latent_representation 'latent_GVAE' \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --homolog_file ${homo}

# 3. generate ldscore
# This strategy uses TSS to assign GSS to SNPs.
hm3_SNP="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/hapmap3_snps/hm"
hg37_annot="/public/share/wchirdzhq2022/Wulab_share/gsMap/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
bfile="/public/share/wchirdzhq2022/Wulab_share/gsMap/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ST_name="E16.5_E2S11.MOSTA"
workdir="/public/share/wchirdzhq2022/Wulab_share/gsMap/Mouse_Embryo"
cmd="gsmap run_generate_ldscore \
        --workdir ${workdir} \
        --sample_name ${ST_name} \
        --chrom {TASK_ID} \
        --bfile_root ${bfile} \
        --keep_snp_root ${hm3_SNP} \
        --gtf_annotation_file ${hg37_annot} \
        --gene_window_size 50000"
qsubshcom "$cmd" 10 100G gsmap 90:00:00 "-array=1-22"

# 4: spatial LDSC
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
hm3_w="/public/share/wchirdzhq2022/Wulab_share/gsMap/LDSC_resource/weights_hm3_no_hla/weights."
cmd="gsmap run_spatial_ldsc \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --sumstats_file ${gwas} \
    --w_file ${hm3_w} \
    --num_processes 30 > spatial_ldsc_condition.log 2>&1"
qsubshcom "$cmd" 10 100G gsmap 10:00:00 ""


# 5: cauchy combination
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
cmd="gsmap run_cauchy_combination \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --annotation 'annotation'"
qsubshcom "$cmd" 1 100G gsmap 10:00:00 ""

# 6: report result
workdir="/public/home/shilulu/Wulab/gsMap/Mouse_Embryo/output"
ST_name="E16.5_E2S11.MOSTA"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/gsMap/ARHL_MVP_AJHG_BBJ_reformatMETAL.sumstats.gz"
cmd="gsmap run_report \
    --workdir ${workdir} \
    --sample_name ${ST_name} \
    --trait_name 'ARHL' \
    --annotation 'annotation' \
    --sumstats_file ${gwas} \
    --top_corr_genes 50"
qsubshcom "$cmd" 1 100G gsmap 10:00:00 ""