#decompress all files.

#modidy file name from Mel. 
for f in A6*.fq; do mv "$f" "${f//A6/Am1}"; done
for f in A8*.fq; do mv "$f" "${f//A8/Am2}"; done
for f in A11*.fq; do mv "$f" "${f//A11/Am5}"; done
for f in A12*.fq; do mv "$f" "${f//A12/Am4}"; done
for f in A17*.fq; do mv "$f" "${f//A17/Am6}"; done

#modify script of Trinity.

## total number of trimmed paired end reads: forward1(unzipped): 3,113,800,000 reads, reverse2(unzipped): 3,122,540,000 reads
## normilisation is a good idea. 

module add UHTS/Assembler/trinityrnaseq/2.4.0
#Trinity version: v2.4.0

Trinity --seqType fq --max_memory 300G \
--left  Am1_434_L5_pairedR1.fastq,Am2_231_L1_pairedR1.fastq,Am2_233_L1_pairedR1.fastq,Am2_274_L1_pairedR1.fastq,Am2_275_L1_pairedR1.fastq,Am2_314_L3_pairedR1.fastq,Am2_314_L5_pairedR1.fastq,Am2_315_L6_pairedR1.fastq,Am2_316_L3_pairedR1.fastq,Am2_318_L2_pairedR1.fastq,Am2_433_L3_pairedR1.fastq,Am2_433_L5_pairedR1.fastq,Am2_434_L4_pairedR1.fastq,Am2_463_L3_pairedR1.fastq,Am2_464_L2_pairedR1.fastq,Am4_231_L3_pairedR1.fastq,Am4_232_L4_pairedR1.fastq,Am4_271_L2_pairedR1.fastq,Am4_272_L3_pairedR1.fastq,Am4_314_L7_pairedR1.fastq,Am4_317_L7_pairedR1.fastq,Am4_434_L3_pairedR1.fastq,Am4_435_L4_pairedR1.fastq,Am4_461_L7_pairedR1.fastq,Am5_231_pairedR1.fastq,Am5_232_L6_pairedR1.fastq,Am5_273_L4_pairedR1.fastq,Am5_274_L5_pairedR1.fastq,Am5_311_L6_pairedR1.fastq,Am5_312_L8_pairedR1.fastq,Am5_313_L7_pairedR1.fastq,Am5_433_L8_pairedR1.fastq,Am5_461_L7_pairedR1.fastq,Am6_462_L6_pairedR1.fastq,Am6_464_L8_pairedR1.fastq \
--right Am1_434_L5_pairedR2.fastq,Am2_231_L1_pairedR2.fastq,Am2_233_L1_pairedR2.fastq,Am2_274_L1_pairedR2.fastq,Am2_275_L1_pairedR2.fastq,Am2_314_L3_pairedR2.fastq,Am2_314_L5_pairedR2.fastq,Am2_315_L6_pairedR2.fastq,Am2_316_L3_pairedR2.fastq,Am2_318_L2_pairedR2.fastq,Am2_433_L3_pairedR2.fastq,Am2_433_L5_pairedR2.fastq,Am2_434_L4_pairedR2.fastq,Am2_463_L3_pairedR2.fastq,Am2_464_L2_pairedR2.fastq,Am4_231_L3_pairedR2.fastq,Am4_232_L4_pairedR2.fastq,Am4_271_L2_pairedR2.fastq,Am4_272_L3_pairedR2.fastq,Am4_314_L7_pairedR2.fastq,Am4_317_L7_pairedR2.fastq,Am4_434_L3_pairedR2.fastq,Am4_435_L4_pairedR2.fastq,Am4_461_L7_pairedR2.fastq,Am5_231_pairedR2.fastq,Am5_232_L6_pairedR2.fastq,Am5_273_L4_pairedR2.fastq,Am5_274_L5_pairedR2.fastq,Am5_311_L6_pairedR2.fastq,Am5_312_L8_pairedR2.fastq,Am5_313_L7_pairedR2.fastq,Am5_433_L8_pairedR2.fastq,Am5_461_L7_pairedR2.fastq,Am6_462_L6_pairedR2.fastq,Am6_464_L8_pairedR2.fastq \
--SS_lib_type RF --normalize_reads --CPU 30 --min_kmer_cov 2 --output Trinity_out_Rtemp 2>&1 | tee run.log