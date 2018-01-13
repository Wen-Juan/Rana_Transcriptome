module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_adults_ \
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



A12FL1_pairedR1.fastq.gz
A12ML1_pairedR1.fastq.gz
A15ML1_pairedR1.fastq.gz
A17FL1_pairedR1.fastq.gz
A17ML1_pairedR1.fastq.gz

_pairedR1.fastq.gz
A2FL1_pairedR1.fastq.gz
_pairedR1.fastq.gz
_pairedR1.fastq.gz
_pairedR1.fastq.gz
A7FL1_pairedR1.fastq.gz
A8FL1_pairedR1.fastq.gz
_pairedR1.fastq.gz

A8ML1_pairedR1.fastq.gz
_pairedR1.fastq.gz


pairedR1.fastq.gz