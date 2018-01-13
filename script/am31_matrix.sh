module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_g31_ \
--name_sample_by_basedir \
Am2_314_all/abundance.tsv \
Am2_315_L6/abundance.tsv \
Am2_316_L3/abundance.tsv \
Am2_318_L2/abundance.tsv \
Am4_314_L7/abundance.tsv \
Am4_317_L7/abundance.tsv \
Am5_311_L6/abundance.tsv \
Am5_312_L8/abundance.tsv \
Am5_313_L7/abundance.tsv
