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

for f in Am*_pairedR1.fastq
do
kallisto quant -i /scratch/beegfs/weekly/wjma/Amm_RNAseq/rt_amm_all.idx -o ./${f%%_pairedR1.fastq} -b 100 ${f%%_pairedR1.fastq}_pairedR1.fastq ${f%%_pairedR1.fastq}_pairedR2.fastq
done

