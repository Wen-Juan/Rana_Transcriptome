#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o logfile_amm1.out
#BUSB -e logfile_amm1.err
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=50000]"
#BSUB -M 50000000
#BSUB -J kallistoamm1.sh

#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0

for f in A10*_R1_pair.fastq
do
kallisto quant -i /scratch/beegfs/weekly/wjma/Amm_RNAseq/rt_amm_all.idx -o ./${f%%_R1_pair.fastq} -b 100 ${f%%_R1_pair.fastq}_R1_pair.fastq ${f%%_R1_pair.fastq}_R2_pair.fastq
done

