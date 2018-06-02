## Rana_Transcriptome
This projects intends to assemble the most comprehensive transcriptome of Rana temporaria, collecting samples from five developmental stages and three adult tissues. This allows the study of sex-biased gene expression throughout development, and its effect on the rate of gene evolution while accounting for pleiotropic expression, which is known to negatively correlate with evolutionary rate.

## A brief pipeline for this project is as follows:

### Step1: Transcriptome is first assembled with Trinity v2.3.0. Following a few filtering steps, transcripts were filtered with possible transposable element insertions by masking the transcriptome assembly using a custom repeat library for Rana temporaria using RepeatMasker, only retaining transcripts that were at least 75% unmasked.

### Step2: Quantifying abundances of transcripts from RNA-Seq with Kallisto v0.43.0.

### Step 3 Analyse differential (sex-biased) gene expression with edgeR v3.4.

### Step 4 Gene divergence analysis (dN/dS), with PRANK v140603 and codeml in PAML.

### Step 5 Tissue specificity Tau calculation, investigating relationship between Tau with dN/dS, and between Tau and sex bias.
