module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_gonad_ \
--name_sample_by_basedir \
A10FO1/abundance.tsv \
A15MT1/abundance.tsv \
A17FO1/abundance.tsv \
A17MT1/abundance.tsv \
A2FO2/abundance.tsv \
A6FO1/abundance.tsv \
A6MT1/abundance.tsv \
A8FO2/abundance.tsv \
A8MT1/abundance.tsv \
A12MT1/abundance.tsv
