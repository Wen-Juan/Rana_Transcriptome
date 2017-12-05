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
--name_sample_by_basedir \
Am1_434_L5/abundance.tsv \
Am2_231_L1/abundance.tsv \
Am2_233_L1/abundance.tsv \
Am2_274_L1/abundance.tsv \
Am2_275_L1/abundance.tsv \
Am2_314_L5/abundance.tsv \
Am2_315_L6/abundance.tsv \
Am2_316_L3/abundance.tsv \
Am2_318_L2/abundance.tsv \
Am2_433_L5/abundance.tsv \
Am2_434_L4/abundance.tsv \
Am2_463_L3/abundance.tsv \
Am2_464_L2/abundance.tsv \
Am4_231_L3/abundance.tsv \
Am4_232_L4/abundance.tsv \
Am4_271_L2/abundance.tsv \
Am4_272_L3/abundance.tsv \
Am4_314_L7/abundance.tsv \
Am4_317_L7/abundance.tsv \
Am4_434_L3/abundance.tsv \
Am4_435_L4/abundance.tsv \
Am4_461_L7/abundance.tsv \
Am5_231/abundance.tsv \
Am5_232_L6/abundance.tsv \
Am5_273_L4/abundance.tsv \
Am5_274_L5/abundance.tsv \
Am5_311_L6/abundance.tsv \
Am5_312_L8/abundance.tsv \
Am5_313_L7/abundance.tsv \
Am5_433_L8/abundance.tsv \
Am5_461_L7/abundance.tsv \
Am6_462_L6/abundance.tsv \
Am6_464_L8/abundance.tsv \
A6MT1/abundance.tsv \
A6FO1/abundance.tsv \
A2FO2/abundance.tsv \
A2FL1/abundance.tsv \
A17MT1/abundance.tsv \
A17ML1/abundance.tsv \
A17MB/abundance.tsv \
A17FO1/abundance.tsv \
A17FL1/abundance.tsv \
A17FB/abundance.tsv \
A16MB/abundance.tsv \
A16FB/abundance.tsv \
A15MT1/abundance.tsv \
A15ML1/abundance.tsv \
A15MB/abundance.tsv \
A15FB/abundance.tsv \
A12MT1/abundance.tsv \
A12ML1/abundance.tsv \
A12FL1/abundance.tsv \
A12FB/abundance.tsv \
A10MB/abundance.tsv \
A10FO1/abundance.tsv \
A10FB/abundance.tsv \
A7FL1/abundance.tsv \
A8FL1/abundance.tsv \
A8FO2/abundance.tsv \
A8MB/abundance.tsv \
A8ML1/abundance.tsv \
A8MT1/abundance.tsv

Step3: #mapping one XX and one XY individuals with highest reads, and then covert SNP positions to the reference transcriptome to generate one XX and one XY transcriptomes.
#mapping with bwa aln
module add UHTS/Aligner/bwa/0.7.13
bwa aln -n 0.04 -o 1 -d 12 -t 20 -M 3 AmmMT_HIMin1.fasta Am4_461_L7_pairedR1.fastq > ./Am4_461_L7_pairedR1.sai
bwa aln -n 0.04 -o 1 -d 12 -t 20 -M 3 AmmMT_HIMin1.fasta Am4_461_L7_pairedR2.fastq > ./Am4_461_L7_pairedR2.sai

module add UHTS/Aligner/bwa/0.7.13
bwa sampe -a 500 -n 3 -N 10 -P AmmMT_HIMin1.fasta Am4_461_L7_pairedR1.sai Am4_461_L7_pairedR2.sai Am4_461_L7_pairedR1.fastq Am4_461_L7_pairedR2.fastq > ./Am4_461_L7.sam

module add UHTS/Analysis/samtools/1.3
samtools view -q 20 -bS Am4_461_L7.sam | samtools sort -o Am4_461_L7.bam

#conver SNPs to the positions in the reference transcriptome
module add UHTS/Analysis/samtools/1.3
samtools mpileup -d8000 -q 25 -Q 30 -uf AmmMT_HIMin1.fasta Am4_461_L7.bam | bcftools call -c | vcfutils.pl vcf2fq > Amm_transcriptome_XY_v2.fq
