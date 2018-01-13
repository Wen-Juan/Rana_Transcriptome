module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_g27_ \
--name_sample_by_basedir \
Am2_274_L1/abundance.tsv \
Am2_275_L1/abundance.tsv \
Am4_271_L2/abundance.tsv \
Am4_272_L3/abundance.tsv \
Am5_273_L4/abundance.tsv \
Am5_274_L5/abundance.tsv
