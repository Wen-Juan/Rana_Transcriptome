#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o ammmap_output.txt
#BUSB -e ammmap_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000000
#BSUB -J bowtie2_map2.sh

module add UHTS/Aligner/bowtie2/2.3.0

bowtie2 -x Amm_devonly_index -q -1 left.fq -2 right.fq -S ammdev.sam 1> amm_devonly.sam 2> bowtie_ammdevonly.log 
