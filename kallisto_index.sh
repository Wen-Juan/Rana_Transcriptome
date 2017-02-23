#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o index_output.txt
#BSUB -e index_error.txt
#BSUB -u wenjuan.ma@unil.ch
#BSUB -N
#BSUB -n 10
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=20000]"
#BSUB -M 20000000
#BSUB -J index.sh

#loading necessary modules and databases
module add UHTS/Analysis/kallisto/0.43.0

kallisto index -i rt_amm_devonly.idx /scratch/beegfs/monthly/wjma/rtRNAseq/all/ref_tseries/Amm_devonly.fasta
