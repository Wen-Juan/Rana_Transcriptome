module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_g46_ \
--name_sample_by_basedir \
Am2_463_L3/abundance.tsv \
Am2_464_L2/abundance.tsv \
Am4_461_L7/abundance.tsv \
Am5_461_L7/abundance.tsv \
Am6_462_L6/abundance.tsv \
Am6_464_L8/abundance.tsv
