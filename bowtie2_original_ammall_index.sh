#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o ammlogfile.out
#BUSB -e ammlogfile.err
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000000
#BSUB -J bowtie2.sh

module add UHTS/Aligner/bowtie2/2.3.0

bowtie2-build Amm_devonly.fasta Amm_devonly_index
