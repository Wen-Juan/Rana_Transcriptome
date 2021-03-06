### venns in R
install.packages("VennDiagram")

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library("SuperExactTest")

#for the shared female-biased genes across five stages.
setwd <-'~/Rana_Transcriptome/'
datapath <- '~/Rana_Transcriptome/input/'
results <- '~/Rana_Transcriptome/output/'

#for tissues from developmental stages, shared female-biased genes
all_sbias <- read.table("~/Rana_Transcriptome/input/venn_diag/de_0.05_1_Am_alltissues.txt", header = TRUE)
str(all_sbias)

all_ortholog_sbias <- read.table("~/Rana_Transcriptome/input/venn_diag/de_0.05_1_orthol_sbun_fi.txt", header = TRUE)
str(all_ortholog_sbias)

###orhologs based on venn diagram for sex-biased gene expression
sbias_g23 <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='G23']
sbias_g27 <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='G27']
sbias_g31 <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='G31']
sbias_g43 <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='G43']
sbias_g46 <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='G46']

venn.plot0 <- venn.diagram(list(G23 = as.character(sbias_g23), G27 = as.character(sbias_g27), G31 = as.character(sbias_g31), G43=as.character(sbias_g43), G46=as.character(sbias_g46)), filename =NULL,
                           fill=c("orange","red","black","green","blue"),
                           ext.line.lwd = 3,
                           cex = 2,
                           cat.cex = 2,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot0),ncol = 1 )

venn_sbias_out0 <- arrangeGrob(gTree(children=venn.plot0),ncol = 1 )
ggsave(file="shared_sbias_larva_orthologs.pdf", venn_sbias_out0, path = "~/Rana_Transcriptome/output/figures/")


sbias_brain <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='Brain']
sbias_gonad <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='Gonad']
sbias_liver <-all_ortholog_sbias$gid[all_ortholog_sbias$stage=='Liver']

venn.plot1 <- venn.diagram(list(Gonad = as.character(sbias_gonad), Brain = as.character(sbias_brain), Liver = as.character(sbias_liver)), filename =NULL,
                           fill=c("orange","green","blue"),
                           ext.line.lwd = 3,
                           cex = 2,
                           cat.cex = 2,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot1),ncol = 1 )

venn_sbias_out1 <- arrangeGrob(gTree(children=venn.plot1),ncol = 1 )
ggsave(file="shared_sbias_adult_orthologs.pdf", venn_sbias_out1, path = "~/Rana_Transcriptome/output/figures/")

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
ggsave(file="shared_sbias_larva.pdf", venn_fbias_out0, path = "~/Rana_Transcriptome/output/figures/")

##investigate shared sex-biased genes between G43 and G46 and the direction of sex bias.
sbias_g43 <- larva_sbias$gid[larva_sbias$stage=='G43']
sbias_g46 <- larva_sbias$gid[larva_sbias$stage=='G46']

venn.plot1 <- venn.diagram(list(G43=as.character(sbias_g43), G46=as.character(sbias_g46)), filename =NULL,
                           fill=c("green","blue"),
                           ext.line.lwd = 3,
                           cex = 2,
                           cat.cex = 2,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot1),ncol = 1 )
venn_fbias_out0 <- arrangeGrob(gTree(children=venn.plot1),ncol = 1 )
ggsave(file="shared_sbias_g43g46.pdf", venn_fbias_out0, path = "~/Rana_Transcriptome/output/figures/")

list_g43g46 <- list(G43=as.character(sbias_g43), G46=as.character(sbias_g46))
overlap_g43g46 <- calculate.overlap(list_g43g46)
overlap_g43g46$a3


##
#test whether the overlap is greater than by chance.
share_5stage <- list(as.character(sbias_g23), as.character(sbias_g27), as.character(sbias_g31), as.character(sbias_g43), as.character(sbias_g46))
list(share_5stage)
str(share_5stage)

total3 <- 4923

####### TO DO for all interections

res=supertest(share_5stage, n=total3)
plot(res, Layout="landscape", degree=2:5, sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Rana_Transcriptome/output/figures/5stage_sb_testbychance.csv", row.names=FALSE)
##
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
ggsave(file="shared_fbias_larva.pdf", venn_fbias_out, path = "~/Rana_Transcriptome/output/figures/")


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
ggsave(file="shared_mbias_larva.pdf", venn_fbias_out1, path = "~/Rana_Transcriptome/output/figures/")

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
ggsave(file="shared_sbias_adult.pdf", venn_fbias_out2, path = "~/Rana_Transcriptome/output/figures/")

#extract overlap list between gonad and liver tissues
gonad_liver <- list(Gonad = as.character(adult_gonad), Liver = as.character(adult_liver))
goliver_overl <- calculate.overlap(gonad_liver)
goliver_overl$a3

gonad_brain <- list(Gonad = as.character(adult_gonad), Brain = as.character(adult_brain))
gobrain_overl <- calculate.overlap(gonad_brain)
gobrain_overl$a3

gonad_G46 <- list(Gonad = as.character(adult_gonad), G46=as.character(sbias_g46))
gobrain_G46 <- calculate.overlap(gonad_G46)
gonad_G46_list <- gobrain_G46$a3

write.table(gonad_G46_list, file="~/Rana_Transcriptome/input/venn_diag/gonadg46_sb_overlap.txt", row.names=FALSE)

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
ggsave(file="shared_fbias_adult.pdf", venn_fbias_out3, path = "~/Rana_Transcriptome/output/figures/")

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
ggsave(file="shared_mbias_adult.pdf", venn_fbias_out4, path = "~/Rana_Transcriptome/output/figures/")


#################
################# G46 vs gonad (the two tissues with high number of sex-biased genes)
#################

################
#######G46 and gonad shared sex-biased genes
################
adult_sbias <- subset(all_sbias,all_sbias$tissue!='larva')
adult_sgonad <- adult_sbias$gid[adult_sbias$stage=='Gonad']
str(adult_sgonad)
larva_sbias <- subset(all_sbias,all_sbias$tissue!='adult')
sbias_g46 <- larva_sbias$gid[larva_sbias$stage=='G46']
head(sbias_g46)

venn.plot5 <- venn.diagram(list(Gonad = as.character(adult_sgonad), G46 = as.character(sbias_g46)), filename =NULL,
                           fill=c("green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5)

grid.arrange(gTree(children=venn.plot5),ncol = 1 )

venn_fbias_out5 <- arrangeGrob(gTree(children=venn.plot5),ncol = 1 )
ggsave(file="shared_sbias_G46gonad.pdf", venn_fbias_out5, path = "~/Rana_Transcriptome/output/figures/")

################
#######G46 and gonad shared female-biased genes
################
adult_fbias <- subset(adult_sbias,adult_sbias$sexbias!='male')
adult_fgonad <- adult_fbias$gid[adult_fbias$stage=='Gonad']
larva_fbias <- subset(larva_sbias,larva_sbias$sexbias!='male')
fbias_g46 <- larva_fbias$gid[larva_fbias$stage=='G46']

venn.plot6 <- venn.diagram(list(Gonad = as.character(adult_fgonad), G46 = as.character(fbias_g46)), filename =NULL,
                           fill=c("green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5)

grid.arrange(gTree(children=venn.plot6),ncol = 1 )

venn_fbias_out6 <- arrangeGrob(gTree(children=venn.plot6),ncol = 1 )
ggsave(file="shared_fbias_G46gonad.pdf", venn_fbias_out6, path = "~/Rana_Transcriptome/output/figures/")

################
#######G46 and gonad shared male-biased genes
################

venn.plot7 <- venn.diagram(list(Gonad = as.character(adult_mgonad), G46 = as.character(mbias_g46)), filename =NULL,
                           fill=c("green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5)

grid.arrange(gTree(children=venn.plot7),ncol = 1 )

venn_fbias_out7 <- arrangeGrob(gTree(children=venn.plot7),ncol = 1 )
ggsave(file="shared_mbias_G46gonad.pdf", venn_fbias_out7, path = "~/Rana_Transcriptome/output/figures/")


#################
################# shared sex-biased genes between G43, G46 and adult tissues.
#################
venn.plot8 <- venn.diagram(list(G43 = as.character(sbias_g43), G46 = as.character(sbias_g46), Gonad = as.character(adult_gonad), Liver=as.character(adult_liver), Brain=as.character(adult_brain)), filename =NULL,
                           fill=c("orange","red","black","green","blue"),
                           ext.line.lwd = 3,
                           cex = 1.5,
                           cat.cex = 1.5,
                           rotation.degree = 65)

grid.arrange(gTree(children=venn.plot8),ncol = 1 )

venn_fbias_out8 <- arrangeGrob(gTree(children=venn.plot8),ncol = 1 )
ggsave(file="shared_sbias_g43g46adults.pdf", venn_fbias_out8, path = "~/Rana_Transcriptome/output/figures/")


#test whether the overlap is greater than by chance.
share_gonadG46 <- list(as.character(adult_gonad), as.character(sbias_g46))
list(share_gonadG46)
str(share_gonadG46)

total3 <- 15591

####### TO DO for all interections

res=supertest(share_gonadG46, n=total3)
plot(res, Layout="landscape", sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Rana_Transcriptome/output/figures/gonadg46_sb_testbychance.csv", row.names=FALSE)

#correction for multiple test.
Pvals = c(0.000013, 0.0002, 0.04,0.4) ### vector of pvals
p.adjust(Pvals, method = "BH") 





