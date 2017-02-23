#1. 
###########assembly statistics###############
794247 transcripts in total

#####Examine assembly stats (not very useful, more for historical reasons)
/software/UHTS/Assembler/trinityrnaseq/2.1.1/util/TrinityStats.pl Amm_devonly.fasta

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  524541
Total trinity transcripts:      794247
Percent GC: 43.61

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 4081
        Contig N20: 2730
        Contig N30: 1966
        Contig N40: 1434
        Contig N50: 1039

        Median contig length: 347
        Average contig: 652.30
        Total assembled bases: 518089483
#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################
        Contig N10: 3289
        Contig N20: 2019
        Contig N30: 1350
        Contig N40: 942
        Contig N50: 675

        Median contig length: 312
        Average contig: 526.85
        Total assembled bases: 276353292 

###################################################
#2 assessing the read contenet of the transcriptome assembly using mapping approach bowtie2.
module add UHTS/Aligner/bowtie2/2.3.0
bowtie2-build Amm_devonly.fasta Amm_devonly.fasta
bowtie2 --local --no-unal -x Amm_devonly.fasta -q -1 left_reads.fq -2 right_reads.fq 


#3. basic filtering of transcriptome. 

