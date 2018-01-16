#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o quan_outputa6.txt
#BSUB -e quan_errora6.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10000000
#BSUB -J quan1_a6.sh

#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0

for f in A6*_pairedR1.fastq.gz
do
kallisto quant -i amm_88perc.idx -t 10 -o ./output/${f%%_pairedR1.fastq.gz} -b 1000 ${f%%_pairedR1.fastq.gz}_pairedR1.fastq.gz ${f%%_pairedR1.fastq.gz}_pairedR2.fastq.gz
done
