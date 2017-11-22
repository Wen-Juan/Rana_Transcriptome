# Rana_Transcriptome
This projects intends to make a complete transcriptome, including 5 developmental stages and 3 adult tissues. We are focusing on sex-biased genes temporal dynamics through development, tissue-specific sex-biased gene expression, detecting dosage compensation, sex chromosome divergence etc.

# A brief pipeline for this project is as follows:

Step1: Transcriptome assembly # by Mellisa Toups
AmmMT_HIMin1.fasta - full transciptome-masked using Dan db at 75% - 61,835 transcripts
#sample used to assembly:
Froglets:
Am2_231_L1
Am2_274_L1
Am2_314_L3
Am2_433_L3
Am2_463_L3

Male - just 1 male - A8:
A8MB - brain
A8ML - liver
A8MT - testes

Female - two females - A8 (couple to male A8), A12
A12FB - brain
A8FL1 - liver
A8FO2 - ovaries

Step2: Kalliso to quantify gene expression.
#separate developmental stages and adult tissues
Total samples:
developmental stages: 5
stage female male: at least 3 per sex per stage, 35 samples in total.

adults:
9 pairs, 3 tissues, 29 libraries (almost 5 replicates per tissue per sex)

Step2.1: developmental stages

# index the transcriptome
module add UHTS/Analysis/kallisto/0.43.0
kallisto index -i am_transcripts.idx /scratch/beegfs/monthly/wjma/amm_rna/expression/AmmMT_HIMin1.fasta
#count the counts of RNAseq

#for each stage, do the quantification using the following scripts.
module add UHTS/Analysis/kallisto/0.43.0
for f in Am2_23*_pairedR1.fastq
do
kallisto quant -i /scratch/beegfs/monthly/wjma/amm_rna/expression/am_transcripts.idx -o /scratch/beegfs/monthly/wjma/amm_rna/expression/kallisto_am/${f%%_pairedR1.fastq} -b 1000 ${f%%_pairedR1.fastq}_pairedR1.fastq ${f%%_pairedR1.fastq}_pairedR2.fastq
done

...
Step2.2: #adult tissues for quantifying gene expression.

Step3: #gene expression matrix generation.
module add UHTS/Assembler/trinityrnaseq/2.4.0
module add UHTS/Analysis/kallisto/0.43.0
module add R/3.3.2

/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto  --out_prefix Amm_ \
Am2_231_L1/abundance.tsv \
Am2_233_L1/abundance.tsv \
Am2_274_L1/abundance.tsv \
...
