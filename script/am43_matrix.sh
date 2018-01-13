module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_g43_ \
--name_sample_by_basedir \
Am1_434_L5/abundance.tsv \
Am2_433_all/abundance.tsv \
Am2_434_L4/abundance.tsv \
Am4_434_L3/abundance.tsv \
Am4_435_L4/abundance.tsv \
Am5_433_L8/abundance.tsv
