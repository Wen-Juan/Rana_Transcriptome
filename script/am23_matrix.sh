module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_g23_ \
--name_sample_by_basedir \
Am2_231_L1/abundance.tsv \
Am2_233_L1/abundance.tsv \
Am4_231_L3/abundance.tsv \
Am4_232_L4/abundance.tsv \
Am5_231/abundance.tsv \
Am5_232_L6/abundance.tsv
