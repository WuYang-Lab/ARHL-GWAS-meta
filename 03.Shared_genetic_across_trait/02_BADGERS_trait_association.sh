#====================================================================#
#                            * BADGERS *                             #
#                          * likely PWAS *                           #
# gets the association between UK_biobank traits and complex disease #
#====================================================================#
conda activate mtag

Badgers="/public/home/shilulu/software/BADGERS/BADGERS/BADGERS.py"
weightdb="/public/home/shilulu/software/BADGERS/BADGERS/UK_Biobank_Round2/weight_db"
covariance="/public/home/shilulu/software/BADGERS/BADGERS/UK_Biobank_Round2/cov"
dbname="/public/home/shilulu/software/BADGERS/BADGERS/1738_traits.csv"
# your gwas and out file
gwas="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/Demeta_raw_hm3_Meta_rescale.txt"
out="/public/home/shilulu/Wulab_project/ARHL/NC_sup_test/06.trait_association/badgers_demeta_out.csv"

python $Badgers \
--model_db_path $weightdb \
--covariance $covariance \
--gwas_path $gwas \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--pvalue_column p \
--beta_column b \
--se_column se \
--output_file $out \
--db_name $dbname
