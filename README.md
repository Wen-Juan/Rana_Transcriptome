## Rana_Transcriptome
This projects intends to assemble the most comprehensive transcriptome of Rana temporaria, collecting samples from five developmental stages and three adult tissues. This allows the study of sex-biased gene expression throughout development, and its effect on the rate of gene evolution while accounting for pleiotropic expression, which is known to negatively correlate with evolutionary rate.

## A brief pipeline for this project is as follows:

### Step1: Transcriptome is first assembled with Trinity v2.3.0. Following a few filtering steps, transcripts were filtered with possible transposable element insertions by masking the transcriptome assembly using a custom repeat library for Rana temporaria using RepeatMasker, only retaining transcripts that were at least 75% unmasked.
the custom repeat library can be found in folder /TE_library
TE_repeats.fasta

### Step2: Quantifying abundances of transcripts from RNA-Seq with Kallisto v0.43.0.
see the scripts in folder /script/
amm_kallisto_pipeline.sh
am23_matrix.sh
am27_matrix.sh
am31_matrix.sh
am43_matrix.sh
am46_matrix.sh
amm_goand_matrix.sh
amm_brain_matrix.sh
amm_liver_matrix.sh

### Step 3 Analyse differential (sex-biased) gene expression with edgeR v3.4.
see the scripts in folder /script/
edgeR_main.r.
stacked_barchart.R

### Step 4 Gene divergence analysis (dN/dS), with PRANK v140603 and codeml in PAML.
see the scripts in folder /script/
Amm_dn_ds.R
Prank_and_PAML.sh
for_loop_codeml.sh

### Step 5 Tissue specificity Tau calculation, investigating relationship between Tau with dN/dS, and between Tau and sex bias.
see the scripts in folder /script/
tau.R
5groups_venn.R
Expression_slidingwindow.r

### Step6 Gene ontology analysis
see the scripts in folder /script/
GO_analysis_05.sh
GO_BP.R
GO_CC.R
GO_MF.R

### Step7 Sex bias gene expression in RPKM and related analysis
RPKM.R
