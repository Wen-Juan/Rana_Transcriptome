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

###############4. Evaluate transcriptome quality
1) N50 to get an idea the number of transcript and genes
###filtering transcriptome
most expressed transcript  --  expressed in at least half of total organs and at least half of samples for each sex. 

evaluate final transcriptome.
1) N50
2) Bowtie2 mapping ####important: using 1>file.sam 2>file.log to get a statistic summary of mapping result.
3) BUSCO
4) by expression value using Kallisto.