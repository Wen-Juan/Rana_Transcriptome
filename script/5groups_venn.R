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

#for tissues from developmental stages, shared female-biased genes
all_sbias <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/venn_diag/de_0.05_1_Am_alltissues.txt", header = TRUE)
str(all_sbias)

##############
############## LARVA
##############

#################
#######larva shared sex-biased genes
#################
larva_sbias <- subset(all_sbias,all_sbias$tissue!='adult')
str(larva_sbias)

sbias_g23 <- larva_sbias$gid[larva_sbias$stage=='G23']
sbias_g27 <- larva_sbias$gid[larva_sbias$stage=='G27']
sbias_g31 <- larva_sbias$gid[larva_sbias$stage=='G31']
sbias_g43 <- larva_sbias$gid[larva_sbias$stage=='G43']
sbias_g46 <- larva_sbias$gid[larva_sbias$stage=='G46']

venn.plot0 <- venn.diagram(list(G23 = as.character(sbias_g23), G27 = as.character(sbias_g27), G31 = as.character(sbias_g31), G43=as.character(sbias_g43), G46=as.character(sbias_g46)), filename =NULL,
                           fill=c("orange","red","black","green","blue"),
                           ext.line.lwd = 3,
                           cex = 2,
                           cat.cex = 2,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot0),ncol = 1 )

venn_fbias_out0 <- arrangeGrob(gTree(children=venn.plot0),ncol = 1 )
ggsave(file="shared_sbias_larva.pdf", venn_fbias_out0, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

##################
#######larvae shared female-biased genes
#################
larva_fbias <- subset(larva_sbias,larva_sbias$sexbias!='male')
str(larva_fbias)

fbis_g23 <- larva_fbias$gid[larva_fbias$stage=='G23']
fbis_g27 <- larva_fbias$gid[larva_fbias$stage=='G27']
fbis_g31 <- larva_fbias$gid[larva_fbias$stage=='G31']
fbis_g43 <- larva_fbias$gid[larva_fbias$stage=='G43']
fbis_g46 <- larva_fbias$gid[larva_fbias$stage=='G46']

venn.plot <- venn.diagram(list(G23 = as.character(fbis_g23), G27 = as.character(fbis_g27), G31 = as.character(fbis_g31), G43=as.character(fbis_g43), G46=as.character(fbis_g46)), filename =NULL,
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2,
                          rotation.degree = 65)
## to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

## to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_fbias_larva.pdf", venn_fbias_out, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")


##############
#######larva shared male-biased genes
##############
larva_mbias <- subset(larva_sbias,larva_sbias$sexbias!='female')
str(larva_mbias)

mbias_g23 <- larva_mbias$gid[larva_mbias$stage=='G23']
mbias_g27 <- larva_mbias$gid[larva_mbias$stage=='G27']
mbias_g31 <- larva_mbias$gid[larva_mbias$stage=='G31']
mbias_g43 <- larva_mbias$gid[larva_mbias$stage=='G43']
mbias_g46 <- larva_mbias$gid[larva_mbias$stage=='G46']

venn.plot1 <- venn.diagram(list(G23 = as.character(mbias_g23), G27 = as.character(mbias_g27), G31 = as.character(mbias_g31), G43=as.character(mbias_g43), G46=as.character(mbias_g46)), filename =NULL,
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2,
                          rotation.degree = 65)

grid.arrange(gTree(children=venn.plot1),ncol = 1 )

venn_fbias_out1 <- arrangeGrob(gTree(children=venn.plot1),ncol = 1 )
ggsave(file="shared_mbias_larva.pdf", venn_fbias_out1, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

##############
############## ADULT
##############

#################
#######adult shared sex-biased genes
#################
adult_sbias <- subset(all_sbias,all_sbias$tissue!='larva')
str(adult_sbias)

adult_gonad <- adult_sbias$gid[adult_sbias$stage=='Gonad']
adult_brain <- adult_sbias$gid[adult_sbias$stage=='Brain']
adult_liver <- adult_sbias$gid[adult_sbias$stage=='Liver']


venn.plot2 <- venn.diagram(list(Gonad = as.character(adult_gonad), Brain = as.character(adult_brain), Liver = as.character(adult_liver)), filename =NULL,
                           fill=c("orange","green","blue"),
                           ext.line.lwd = 3,
                           cex = 2,
                           cat.cex = 2,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot2),ncol = 1 )

venn_fbias_out2 <- arrangeGrob(gTree(children=venn.plot2),ncol = 1 )
ggsave(file="shared_sbias_adult.pdf", venn_fbias_out2, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")


#################
#######adult shared female-biased genes
#################
adult_fbias <- subset(adult_sbias,adult_sbias$sexbias!='male')
str(adult_fbias)

adult_fgonad <- adult_fbias$gid[adult_fbias$stage=='Gonad']
adult_fbrain <- adult_fbias$gid[adult_fbias$stage=='Brain']
adult_fliver <- adult_fbias$gid[adult_fbias$stage=='Liver']


venn.plot3 <- venn.diagram(list(Gonad = as.character(adult_fgonad), Brain = as.character(adult_fbrain), Liver = as.character(adult_fliver)), filename =NULL,
                           fill=c("orange","green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot3),ncol = 1 )

venn_fbias_out3 <- arrangeGrob(gTree(children=venn.plot3),ncol = 1 )
ggsave(file="shared_fbias_adult.pdf", venn_fbias_out3, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

#################
#######adult shared male-biased genes
#################
adult_mbias <- subset(adult_sbias,adult_sbias$sexbias!='female')
str(adult_mbias)

adult_mgonad <- adult_mbias$gid[adult_mbias$stage=='Gonad']
adult_mbrain <- adult_mbias$gid[adult_mbias$stage=='Brain']
adult_mliver <- adult_mbias$gid[adult_mbias$stage=='Liver']


venn.plot4 <- venn.diagram(list(Gonad = as.character(adult_mgonad), Brain = as.character(adult_mbrain), Liver = as.character(adult_mliver)), filename =NULL,
                           fill=c("orange","green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot4),ncol = 1 )

venn_fbias_out4 <- arrangeGrob(gTree(children=venn.plot4),ncol = 1 )
ggsave(file="shared_mbias_adult.pdf", venn_fbias_out4, path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")


#################
################# G46 vs gonad (the two tissues with high number of sex-biased genes)
#################

#######adult shared male-biased genes
#################


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





