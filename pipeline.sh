###################1
#decompress all files.

#modidy file name from Mel. 
for f in A6*.fq; do mv "$f" "${f//A6/Am1}"; done
for f in A8*.fq; do mv "$f" "${f//A8/Am2}"; done
for f in A12*.fq; do mv "$f" "${f//A12/Am4}"; done
for f in A17*.fq; do mv "$f" "${f//A17/Am6}"; done
for f in A*.fq; do mv "$f" "${f//fq/fastq}"; done

###############2. trimming raw reads (see pipeline of RNAseq_sexdet)

###############3. de novo assembly using Trinity
#modify script of Trinity.

## total number of trimmed paired end reads: forward1(unzipped): 4,647,130,000 reads, reverse2(unzipped): 4,648,450,000 reads
## normilisation is a good idea. 

module add UHTS/Assembler/trinityrnaseq/2.4.0
#Trinity version: v2.4.0

Trinity --seqType fq --max_memory 300G \
--left  Am1_434_L5_pairedR1.fastq,Am2_231_L1_pairedR1.fastq,Am2_233_L1_pairedR1.fastq,Am2_274_L1_pairedR1.fastq,Am2_275_L1_pairedR1.fastq,Am2_314_L3_pairedR1.fastq,Am2_314_L5_pairedR1.fastq,Am2_315_L6_pairedR1.fastq,Am2_316_L3_pairedR1.fastq,Am2_318_L2_pairedR1.fastq,Am2_433_L3_pairedR1.fastq,Am2_433_L5_pairedR1.fastq,Am2_434_L4_pairedR1.fastq,Am2_463_L3_pairedR1.fastq,Am2_464_L2_pairedR1.fastq,Am4_231_L3_pairedR1.fastq,Am4_232_L4_pairedR1.fastq,Am4_271_L2_pairedR1.fastq,Am4_272_L3_pairedR1.fastq,Am4_314_L7_pairedR1.fastq,Am4_317_L7_pairedR1.fastq,Am4_434_L3_pairedR1.fastq,Am4_435_L4_pairedR1.fastq,Am4_461_L7_pairedR1.fastq,Am5_231_pairedR1.fastq,Am5_232_L6_pairedR1.fastq,Am5_273_L4_pairedR1.fastq,Am5_274_L5_pairedR1.fastq,Am5_311_L6_pairedR1.fastq,Am5_312_L8_pairedR1.fastq,Am5_313_L7_pairedR1.fastq,Am5_433_L8_pairedR1.fastq,Am5_461_L7_pairedR1.fastq,Am6_462_L6_pairedR1.fastq,Am6_464_L8_pairedR1.fastq,A10FB_R1_pair.fastq,A10FO1_R1_pair.fastq,A10MB_R1_pair.fastq,A15FB_R1_pair.fastq,A15MB_R1_pair.fastq,A15ML1_R1_pair.fastq,A15MT1_R1_pair.fastq,A16FB_R1_pair.fastq,A16MB_R1_pair.fastq,A2FL1_R1_pair.fastq,A2FO2_R1_pair.fastq,A7FL1_R1_pair.fastq,Am1FO1_R1_pair.fastq,Am1ML1_R1_pair.fastq,Am1MT1_R1_pair.fastq,Am2FL1_R1_pair.fastq,Am2FO2_R1_pair.fastq,Am2MB_R1_pair.fastq,Am2ML1_R1_pair.fastq,Am2MT1_R1_pair.fastq,Am4FB_R1_pair.fastq,Am4FL1_R1_pair.fastq,Am4ML1_R1_pair.fastq,Am4MT1_R1_pair.fastq,Am6FB_R1_pair.fastq,Am6FL1_R1_pair.fastq,Am6FO1_R1_pair.fastq,Am6MB_R1_pair.fastq,Am6ML1_R1_pair.fastq,Am6MT1_R1_pair.fastq \
--right Am1_434_L5_pairedR2.fastq,Am2_231_L1_pairedR2.fastq,Am2_233_L1_pairedR2.fastq,Am2_274_L1_pairedR2.fastq,Am2_275_L1_pairedR2.fastq,Am2_314_L3_pairedR2.fastq,Am2_314_L5_pairedR2.fastq,Am2_315_L6_pairedR2.fastq,Am2_316_L3_pairedR2.fastq,Am2_318_L2_pairedR2.fastq,Am2_433_L3_pairedR2.fastq,Am2_433_L5_pairedR2.fastq,Am2_434_L4_pairedR2.fastq,Am2_463_L3_pairedR2.fastq,Am2_464_L2_pairedR2.fastq,Am4_231_L3_pairedR2.fastq,Am4_232_L4_pairedR2.fastq,Am4_271_L2_pairedR2.fastq,Am4_272_L3_pairedR2.fastq,Am4_314_L7_pairedR2.fastq,Am4_317_L7_pairedR2.fastq,Am4_434_L3_pairedR2.fastq,Am4_435_L4_pairedR2.fastq,Am4_461_L7_pairedR2.fastq,Am5_231_pairedR2.fastq,Am5_232_L6_pairedR2.fastq,Am5_273_L4_pairedR2.fastq,Am5_274_L5_pairedR2.fastq,Am5_311_L6_pairedR2.fastq,Am5_312_L8_pairedR2.fastq,Am5_313_L7_pairedR2.fastq,Am5_433_L8_pairedR2.fastq,Am5_461_L7_pairedR2.fastq,Am6_462_L6_pairedR2.fastq,Am6_464_L8_pairedR2.fastq,A10FB_R2_pair.fastq,A10FO1_R2_pair.fastq,A10MB_R2_pair.fastq,A15FB_R2_pair.fastq,A15MB_R2_pair.fastq,A15ML1_R2_pair.fastq,A15MT1_R2_pair.fastq,A16FB_R2_pair.fastq,A16MB_R2_pair.fastq,A2FL1_R2_pair.fastq,A2FO2_R2_pair.fastq,A7FL1_R2_pair.fastq,Am1FO1_R2_pair.fastq,Am1ML1_R2_pair.fastq,Am1MT1_R2_pair.fastq,Am2FL1_R2_pair.fastq,Am2FO2_R2_pair.fastq,Am2MB_R2_pair.fastq,Am2ML1_R2_pair.fastq,Am2MT1_R2_pair.fastq,Am4FB_R2_pair.fastq,Am4FL1_R2_pair.fastq,Am4ML1_R2_pair.fastq,Am4MT1_R2_pair.fastq,Am6FB_R2_pair.fastq,Am6FL1_R2_pair.fastq,Am6FO1_R2_pair.fastq,Am6MB_R2_pair.fastq,Am6ML1_R2_pair.fastq,Am6MT1_R2_pair.fastq \
--SS_lib_type RF --normalize_reads --CPU 30 --min_kmer_cov 2 --output Trinity_out_Rtamm_all 2>&1 | tee run.log

###############4. Evaluate transcriptome qualit
###before filtering#########################################
1.1) N50 to get an idea the number of transcript and genes
############################################################
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  571975
Total trinity transcripts:      1078742
Percent GC: 44.30

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3847
        Contig N20: 2499
        Contig N30: 1770
        Contig N40: 1268
        Contig N50: 910

        Median contig length: 351
        Average contig: 624.52
        Total assembled bases: 673696009


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3145
        Contig N20: 1936
        Contig N30: 1296
        Contig N40: 920
        Contig N50: 677

        Median contig length: 330
        Average contig: 538.29
        Total assembled bases: 307888670
######################################
1.2) remove smaller-length transcripts (<300bp)
1.3) reassess the N50 statistics.


2) mapping
3) expression using Kallisto

###filtering transcriptome
most expressed transcript  --  expressed in at least half of total organs and at least half of samples for each sex. 

evaluate final transcriptome.
1) N50
2) Bowtie2 mapping ####important: using 1>file.sam 2>file.log to get a statistic summary of mapping result.
3) BUSCO
4) by expression value using Kallisto.

####################4######################
####mapping reads to TE library generated from XY Spain male#################
####4.1: Sex specific TE? Ammarnas: brain XXfemale vs. XYmale

step1: make an index of the target library
module add UHTS/Analysis/kallisto/0.43.0
kallisto index -i rt_TE_spain.idx /scratch/beegfs/weekly/wjma/Amm_RNAseq/rana_TEs.fasta

step2: kallisto quantify the expression on TEs if there are any.
#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0
###There are 5 female brains and 5 male brains.

for f in A10*_R1_pair.fastq
do
kallisto quant -i /scratch/beegfs/weekly/wjma/Amm_RNAseq/rt_TE_spain.idx -o ./${f%%_R1_pair.fastq} -b 100 ${f%%_R1_pair.fastq}_R1_pair.fastq ${f%%_R1_pair.fastq}_R2_pair.fastq
done

####4.2: Tissue specific in XY males TE? Ammarnas: XY brain vs XY liver vs XY gonad
#similar script as before. 

####4.3: TE different among XY, XY', XX? Amm_XY_froglets vs. Tvedora_XY'_froglet
#Tv1_432_L2.fastq
#Tv3_431_L7.fastq
#Tv2_464_L3.fastq
#Tv1_463_L6.fastq
#Tv2_461_L3.fastq
#Tv2_461_L5.fastq
#Am2_434_L4.fastq
#Am1_434_L5.fastq
#Am4_434_L3.fastq
#Am4_461_L7.fastq
#Am2_463_L3.fastq
#Am6_462_L6.fastq

#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o logfile_ammTv1.out
#BUSB -e logfile_ammTv1.err
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J kallistoAmmTv1.sh

#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0
###There are 6 tv male froglets and 6 ammarnas male froglets

for f in Tv1*_pairedR1.fastq
do
kallisto quant -i /scratch/beegfs/weekly/wjma/Amm_RNAseq/rt_TE_spain.idx -o ./${f%%_pairedR1.fastq} -b 100 ${f%%_pairedR1.fastq}_pairedR1.fastq ${f%%_pairedR1.fastq}_pairedR2.fastq
done

#####5.expression count matrix and differential expression analysis######################
###############comparison between Ammarnas XX female and XY male brains 
module add UHTS/Assembler/trinityrnaseq/2.4.0
module add R/3.3.2
/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --out_prefix Amm_FM_brain \
--name_sample_by_basedir \
A10FB/abundance.tsv \
A10MB/abundance.tsv \
A15FB/abundance.tsv \
A15MB/abundance.tsv \
A16FB/abundance.tsv \
A16MB/abundance.tsv \
Am2MB/abundance.tsv \
Am4FB/abundance.tsv \
Am6FB/abundance.tsv \
Am6MB/abundance.tsv

###############comparison among Ammarnas XY brain, liver and testis
module add UHTS/Assembler/trinityrnaseq/2.4.0
module add R/3.3.2
/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --out_prefix Amm_XY_BLT \
--name_sample_by_basedir \
A10MB/abundance.tsv \
A15MB/abundance.tsv \
A15ML1/abundance.tsv \
A15MT1/abundance.tsv \
A16MB/abundance.tsv \
Am1MT1/abundance.tsv \
Am2MB/abundance.tsv \
Am2ML1/abundance.tsv \
Am2MT1/abundance.tsv \
Am4ML1/abundance.tsv \
Am4MT1/abundance.tsv \
Am6MB/abundance.tsv \
Am6ML1/abundance.tsv \
Am6MT1/abundance.tsv

###############comparison among Ammarnas XY brain, liver and testis
module add UHTS/Assembler/trinityrnaseq/2.4.0
module add R/3.3.2
/software/UHTS/Assembler/trinityrnaseq/2.4.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --out_prefix Amm_TV_froglets \
--name_sample_by_basedir \
Am1_434_L5/abundance.tsv \
Am2_434_L4/abundance.tsv \
Am2_463_L3/abundance.tsv \
Am4_434_L3/abundance.tsv \
Am4_461_L7/abundance.tsv \
Am6_462_L6/abundance.tsv \
Tv1_432_L2/abundance.tsv \
Tv1_463_L6/abundance.tsv \
Tv2_461_L3/abundance.tsv \
Tv2_461_L5/abundance.tsv \
Tv2_464_L3/abundance.tsv \
Tv3_431_L7/abundance.tsv























