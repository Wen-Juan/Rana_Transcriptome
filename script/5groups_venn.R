### venns in R
install.packages("VennDiagram")

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#for the shared female-biased genes across five stages.
setwd <-'/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/'
datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/'
results <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/'

#for shared female-biased genes
fivestages_fbias <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/female_bias_stages.txt", header = TRUE)
str(fivestages_fbias)

adults_fbias <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/female_bias_adults.txt", header = TRUE)
str(adults_fbias)

#fbis_g23 <- fivestages_fbias$trans[fivestages_fbias$stage=='G23']
fbis_brain <- adults_fbias$trans[adults_fbias$stage=='Brain']
#fbis_g27 <- fivestages_fbias$trans[fivestages_fbias$stage=='G27']
fbis_ovary <- adults_fbias$trans[adults_fbias$stage=='Ovary']
#fbis_g31 <- fivestages_fbias$trans[fivestages_fbias$stage=='G31']
fbis_liver <- adults_fbias$trans[adults_fbias$stage=='Liver']
fbis_g43 <- adults_fbias$trans[adults_fbias$stage=='G43']
fbis_g46 <- adults_fbias$trans[adults_fbias$stage=='G46']

venn.plot <- venn.diagram(list(G23 = as.character(fbis_g23), G27 = as.character(fbis_g27), G31 = as.character(fbis_g31), G43=as.character(fbis_g43), G46=as.character(fbis_g46)), filename =NULL,
                          #fill=rainbow(5),
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2,
                          rotation.degree = 65)


venn.plot <- venn.diagram(list(Brain = as.character(fbis_brain), Ovary = as.character(fbis_ovary), Liver = as.character(fbis_liver), G43=as.character(fbis_g43), G46=as.character(fbis_g46)), filename =NULL,
                          #fill=rainbow(5),
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 1.5,
                          rotation.degree = 65)

## to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

## to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_fbias_adults.pdf", venn_fbias_out, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

#for shared male-biased genes
mbias <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/male_bias_stages.txt", header = TRUE)
str(mbias)

mbias_adult <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/male_bias_adults.txt", header = TRUE)
str(mbias_adult)

#mbis_g23 <- mbias$trans[mbias$stage=='G23']
mbis_brain <- mbias_adult$trans[mbias_adult$stage=='Brain']
#mbis_g27 <- mbias$trans[mbias$stage=='G27']
mbis_liver <- mbias_adult$trans[mbias_adult$stage=='Liver']
#mbis_g31 <- mbias$trans[mbias$stage=='G31']
mbis_testis <- mbias_adult$trans[mbias_adult$stage=='Testis']
mbis_g43 <- mbias_adult$trans[mbias_adult$stage=='G43']
mbis_g46 <- mbias_adult$trans[mbias_adult$stage=='G46']

venn.plot.2 <- venn.diagram(list(G23 = as.character(mbis_g23),G27 = as.character(mbis_g27), G31 = as.character(mbis_g31), G43=as.character(mbis_g43), G46=as.character(mbis_g46)), filename =NULL,
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2,
                          rotation.degree = 65)

venn.plot <- venn.diagram(list(Brain = as.character(mbis_brain), Testis = as.character(mbis_testis), Liver = as.character(mbis_liver), G43=as.character(mbis_g43), G46=as.character(mbis_g46)), filename =NULL,
                          #fill=rainbow(5),
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 1.5,
                          rotation.degree = 65)

## to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )
## to output to pdf
venn_mbias_out <- arrangeGrob(gTree(children=venn.plot.2),ncol = 1 )
ggsave(file="shared_mbias_adults.pdf", venn_mbias_out, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

#for shared unbiased genes
fivestages_unbias <- read.table("/Users/Wen-Juan/git/rana_transcriptome/input/sex_bias/share_unbias.txt", header = TRUE)
str(fivestages_unbias)

unbias_g23 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G23']
unbias_g27 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G27']
unbias_g31 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G31']
unbias_g43 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G43']
unbias_g46 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G46']

venn.plot <- venn.diagram(list(G23 = as.character(unbias_g23), G27 = as.character(unbias_g27), G31 = as.character(unbias_g31), G43=as.character(unbias_g43), G46=as.character(unbias_g46)), filename =NULL,
                          cat.col=c("orange","red","black","green","blue"),
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 2.5,
                          rotation.degree = 60)

grid.arrange(gTree(children=venn.plot),ncol = 1 )

venn_unbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_unbias.pdf", venn_unbias_out, path = "/Users/Wen-Juan/git/rana_transcriptome/input/sex_bias/")


#for shared sex-biased and unbiased genes between G43 and G46
g4346_fbias <- read.table("/Users/Wen-Juan/git/rana_transcriptome/input/sex_bias/G4346_unbias_sb.txt", header = TRUE)
str(g4346_fbias)

sb_g43 <- g4346_fbias$trans[g4346_fbias$stage=='G43_sb']
unbis_g43 <- g4346_fbias$trans[g4346_fbias$stage=='G43_unbias']
sb_g46 <- g4346_fbias$trans[g4346_fbias$stage=='G46_sb']
unbis_g46 <- g4346_fbias$trans[g4346_fbias$stage=='G46_unbias']

venn.plot.1 <- venn.diagram(list(G43_SB = as.character(sb_g43),G43_UN = as.character(unbis_g43),G46_SB=as.character(sb_g46),G46_UN=as.character(unbis_g46)), filename =NULL,
                          fill=c("red","green","orange","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 2)

## to draw to the screen:
grid.arrange(gTree(children=venn.plot.1),ncol = 1 )
## to output to pdf
venn_shareg4346_out <- arrangeGrob(gTree(children=venn.plot.1),ncol = 1 )
ggsave(file="G4346_shared.pdf", venn_shareg4346_out, path = "/Users/Wen-Juan/git/rana_transcriptome/input/sex_bias/")





