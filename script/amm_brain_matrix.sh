module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_brain_ \
--name_sample_by_basedir \
A10FB/abundance.tsv \
A10MB/abundance.tsv \
A12FB/abundance.tsv \
A15FB/abundance.tsv \
A15MB/abundance.tsv \
A16FB/abundance.tsv \
A16MB/abundance.tsv \
A17FB/abundance.tsv \
A17MB/abundance.tsv \
A8MB/abundance.tsv
