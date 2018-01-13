module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_gonad_ \
--name_sample_by_basedir \
A12FL1/abundance.tsv \
A12ML1/abundance.tsv \
A15ML1/abundance.tsv \
A17FL1/abundance.tsv \
A17ML1/abundance.tsv \
A2FL1/abundance.tsv \
A7FL1/abundance.tsv \
A8FL1/abundance.tsv \
A8ML1/abundance.tsv
