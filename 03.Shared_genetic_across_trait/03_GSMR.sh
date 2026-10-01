#--------------------------- GSMR ----//
#
# GSMR for some traits select from BADGERS
# 
#--------------------------- GSMR ----//
dir="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS"
echo "Tiredness $dir/clean_tiredness_lethargy2080.fastGWA.gz
Neuroticism $dir/clean_Neuroticism_score.fastGWA.gz
LongStandIllness $dir/clean_Long-standing_illness.fastGWA.gz
Overall_health_rating $dir/clean_Overall_health_rating.fastGWA.gz
Noisy_workplace $dir/clean_Noisy_workplace.fastGWA.gz
Snoring $dir/clean_Snoring.fastGWA.gz
Insomnia $dir/clean_Sleeplessness-insomnia.fastGWA.gz
Taking_other_prescription_medications $dir/clean_Taking_other_prescription_medications.fastGWA.gz
Waist_circumference $dir/Waist_circumference.fastGWA.gz
Loneliness $dir/clean_loneliness.fastGWA.gz
derpess_mood $dir/clean_depressed_mood.fastGWA.gz
Past_tobacco_smoking $dir/clean_Past_tobacco_smoking.fastGWA.gz" > exposure.txt


echo "Waist_circumference $dir/Waist_circumference.fastGWA.gz" > exposure_waist.txt
echo "ARHL /public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/Demeta_raw_hm3_Meta_rescale.txt" > outcome.txt

#=============================================#
#! /bin/bash 
gcta="/public/home/wuyang/bin/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
ref="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
exposure="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR/exposure.txt"
outcome="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/GSMR/outcome.txt"

$gcta --bfile $ref \
--gsmr-file $exposure $outcome \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--diff-freq 0.5 \
--clump-r2 0.05 \
--gsmr2-beta \
--gsmr-snp-min 10 \
--heidi-thresh 0.01 \
--effect-plot \
--out ARHL_gsmr \
--thread-num 20

# no heidi
$gcta --bfile $ref \
--gsmr-file $exposure $outcome \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--diff-freq 0.5 \
--clump-r2 0.05 \
--gsmr2-beta \
--gsmr-snp-min 10 \
--heidi-thresh 1e-700 \
--effect-plot \
--out ARHL_gsmr_noheidi \
--thread-num 20