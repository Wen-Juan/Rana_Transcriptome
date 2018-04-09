# laading R libraries
library(edgeR)
library(gplots)
library(ggplot2)
install.packages("gdata", dependencies=TRUE) #for cbind of multiple data frames
library(gdata)
install.packages("plyr")
library(plyr)
library(ggpubr)

#loading datasets
logFC_23 <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amg23.txt", header = T)
str(logFC_23)
logFC_23$stage = "G23"

logFC_27 <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amg27.txt", header = T)
str(logFC_27)
logFC_27$stage = "G27"

logFC_31 <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amg31.txt", header = T)
str(logFC_31)
logFC_31$stage = "G31"

logFC_43 <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amg43.txt", header = T)
str(logFC_43)
logFC_43$stage = "G43"

logFC_46nosr <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amg46nosr.txt", header = T)
str(logFC_46nosr)
logFC_46nosr$stage = "G46"

logFC_gonad <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amgonad.txt", header = T)
str(logFC_gonad)
logFC_gonad$stage = "gonad"

logFC_liver <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Amliver.txt", header = T)
str(logFC_liver)
logFC_liver$stage = "liver"

logFC_brain <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/LogCPM_0.05_Ambrain.txt", header = T)
str(logFC_brain)
logFC_brain$stage = "brain"

g23_sub <- data.frame(logFC_23$logFC.XY23.XX23,logFC_23$stage)
write.table(g23_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g23_sub.txt", sep="\t", col.names=F)
g27_sub <- data.frame(logFC_27$logFC.XY27.XX27,logFC_27$stage)
write.table(g27_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g27_sub.txt", sep="\t", col.names=F)
g31_sub <- data.frame(logFC_31$logFC.XY31.XX31,logFC_31$stage)
write.table(g31_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g31_sub.txt", sep="\t", col.names=F)
g43_sub <- data.frame(logFC_43$logFC.XY43.XX43,logFC_43$stage)
write.table(g43_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_sub.txt", sep="\t", col.names=F)
g46_sub <- data.frame(logFC_46nosr$logFC.XY46.XX46,logFC_46nosr$stage)
write.table(g46_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_sub.txt", sep="\t", col.names=F)
gonad_sub <- data.frame(logFC_gonad$logFC.XYtestis.XXovary,logFC_gonad$stage)
write.table(gonad_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_sub.txt", sep="\t", col.names=F)
brain_sub <- data.frame(logFC_brain$logFC.XYbrain.XXrain,logFC_brain$stage)
write.table(brain_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/brain_sub.txt", sep="\t", col.names=F)
liver_sub <- data.frame(logFC_liver$logFC.XYliver.XXliver,logFC_liver$stage)
write.table(liver_sub, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_sub.txt", sep="\t", col.names=F)

logFC_all <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/all_sub_logfc.txt", header = T)
str(logFC_all)

new_stage <- factor(logFC_all$stage, levels=c("G23", "G27","G31","G43", "G46", "gonad","liver", "brain"))

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/violin_expressionratio.pdf")
p = ggplot(logFC_all, aes(factor(new_stage),logFC), ylim)
p + geom_violin(scale="width",fill="dodgerblue3", trim=F) + 
  geom_hline(yintercept=-1,linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") +
  theme_set(theme_bw(base_size=12)) +
  #+ geom_hline(yintercept=-2,linetype="dashed") + geom_hline(yintercept=2, linetype="dashed") + 
  #theme(panel.background = element_rect(fill='blue')) + 
  scale_y_continuous(breaks=c(-15,-10, -5, -1, 0, 1, 5, 10, 15)) +theme(axis.text=element_text(size=12), 
#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
dev.off()

geom_boxplot(notch=TRUE)  +
  ylim(-2,5) +
  scale_fill_manual(values = c("red","blue")) +
  theme(axis.title.x = element_text(size=15,colour = "black"),axis.title.y = element_text(size=15,colour = "black")) + 
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

### dataset from G46
datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/rpkm_list_count_sorted.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
kdata1 <- subset(kdata, kdata$logFC.XY46.XX46<=-1 & kdata$FDR.XY46.XX46< 0.05)
kdata2 <- subset(kdata, kdata$logFC.XY46.XX46>=1 & kdata$FDR.XY46.XX46< 0.05)
kdata3 <- subset(kdata, kdata$logFC.XY46.XX46>-1 & kdata$logFC.XY46.XX46<1 & kdata$FDR.XY46.XX46< 0.05)

kdata1$meanXY <- (kdata1$logAm2_463 + kdata1$logAm4_461 + kdata1$logAm6_462) / 3 #XY male group at G46
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$logAm2_463 + kdata2$logAm4_461 + kdata2$logAm6_462) / 3 #XY male group at G46
kdata2$group1="meanmale"

kdata3$meanXY <- (kdata3$logAm2_463 + kdata3$logAm4_461 + kdata3$logAm6_462) / 3 #XY male group at G46
kdata3$group1="meanmale"

kdata1$meanXX <- (kdata1$logAm2_464+kdata1$logAm6_464)/2 #XX female group at G46
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$logAm2_464+kdata2$logAm6_464)/2 #XX female group at G46
kdata2$group2="meanfemale"

kdata3$meanXX <- (kdata3$logAm2_464+kdata3$logAm6_464)/2 #XY female group at G46
kdata3$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46nosr_fbias.txt", sep="\t", col.names=F)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46nosr_mbias.txt", sep="\t", col.names=F)
write.table(kdata3, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46nosr_unbias.txt", sep="\t", col.names=F)

g46_fbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_fbias_all.txt", header=TRUE)
g46_fbias_all$bias = "female"
g46_mbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_mbias_all.txt", header=TRUE)
g46_mbias_all$bias = "male"
g46_unbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_unbias_all.txt", header=TRUE)
g46_unbias_all$bias = "unbias"

write.table(g46_fbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_fbias_all2.txt", sep="\t", col.names=F)
write.table(g46_mbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_mbias_all2.txt", sep="\t", col.names=F)
write.table(g46_unbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_unbias_all2.txt", sep="\t", col.names=F)

g46_bias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_sbunbias_all_fi.txt", header=TRUE)
str(g46_bias_all)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sexbias_logcpmg46_fbias.pdf")

sts <- boxplot.stats(g46_bias_all$log2RPKM)$stats
p2 <- ggplot(data = g46_bias_all, aes(x=Bias, y= log2RPKM, fill=Sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0.8)) +
  coord_cartesian(ylim = c(sts[1]*1.8,max(sts)*1.3)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
  #theme(legend.position="none")
dev.off()

###cut-off figures of G46 Log2FC
datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/rpkm_list_count_sorted.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
kdata1 <- subset(kdata, kdata$logFC.XY46.XX46<=-1 & kdata$FDR.XY46.XX46< 0.05)
kdata2 <- subset(kdata, kdata$logFC.XY46.XX46>=1 & kdata$FDR.XY46.XX46< 0.05)

kdata1$meanXY <- (kdata1$Am2_463 + kdata1$Am4_461 + kdata1$Am6_462) / 3 #XY male group at G46
kdata1$meanXX <- (kdata1$Am2_464 + kdata1$Am6_464)/2

kdata1$log2meanXY <- log2(kdata1$meanXY + 0.0001)
kdata1$group1="meanmale"
kdata1$log2meanXX <- log2(kdata1$meanXX + 0.0001)
kdata1$group1="meanfemale"

kdata2$meanXY <- (kdata2$Am2_463 + kdata2$Am4_461 + kdata2$Am6_462) / 3 #XY male group at G46
kdata2$meanXX <- (kdata2$Am2_464+kdata2$Am6_464)/2 #XX female group at G46

kdata2$log2meanXY <- log2(kdata2$meanXY + 0.0001)
kdata2$group1="meanmale"
kdata2$log2meanXX <- log2(kdata2$meanXX + 0.0001)
kdata2$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_fbias_rpkmfdr.txt", sep="\t", col.names=T)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_mbias_rpkmfdr.txt", sep="\t", col.names=T)

g46_fbias_rpkmfrd_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_fbias_all_cutoff_rpkm_fifi.txt", header = T)
str(g46_fbias_rpkmfrd_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/cutoff_rpkm_g46_fbias.pdf")

p_g46fb <- ggplot(data = g46_fbias_rpkmfrd_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-13,-10,-5,0,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

g46_mbias_rpkmfrd_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_mbias_all_rpkmfdr_cutoff_fi.txt", header = T)
str(g46_mbias_rpkmfrd_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/cutoff_rpkm_g46_mbias.pdf")

p_g46mb <- ggplot(data = g46_mbias_rpkmfrd_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-13,-10,-5,0,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/cutoff_rpkm_g46_fi.pdf")

ggarrange(p_g46fb, p_g46mb, labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

###stats on cut-off sex bias
male_g46_mbias_rpkm <- subset(g46_mbias_rpkmfrd_cutoff,g46_mbias_rpkmfrd_cutoff$sex=='meanmale')
female_g46_mbias_rpkm <- subset(g46_mbias_rpkmfrd_cutoff,g46_mbias_rpkmfrd_cutoff$sex=='meanfemale')

wilcox.test(male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC1'], male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC3'])#W = 12610, p-value = 0.01116
wilcox.test(male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC3'], male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC5'])#W = 48, p-value = 0.8831
wilcox.test(male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC5'], male_g46_mbias_rpkm $log2RPKM[male_g46_mbias_rpkm $cutoff=='FC7'])#W = 29, p-value = 0.8696

wilcox.test(female_g46_mbias_rpkm $log2RPKM[female_g46_mbias_rpkm $cutoff=='FC1'], female_g46_mbias_rpkm $log2RPKM[female_g46_mbias_rpkm $cutoff=='FC3'])#W = 17909, p-value = 6.374e-14
wilcox.test(female_g46_mbias_rpkm $log2RPKM[female_g46_mbias_rpkm $cutoff=='FC3'], female_g46_mbias_rpkm $log2RPKM[female_g46_mbias_rpkm $cutoff=='FC5'])#W = 74, p-value = 0.07478
wilcox.test(female_g46_mbias_rpkm $log2RPKM[female_g46_mbias_rpkm $cutoff=='FC5'], female_g46_mbias_rpkm $log2RPKM[female_g46_fbias_rpkm $cutoff=='FC7']) #W = 0, p-value = 0.003031

male_G46_fbias_rpkm <- subset(g46_fbias_rpkmfrd_cutoff,g46_fbias_rpkmfrd_cutoff$sex=='meanmale')
female_G46_fbias_rpkm <- subset(g46_fbias_rpkmfrd_cutoff,g46_fbias_rpkmfrd_cutoff$sex=='meanfemale')

wilcox.test(male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC1'], male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC3'])#W = 1141500, p-value < 2.2e-16
wilcox.test(male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC3'], male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC5'])#W = 48040, p-value < 2.2e-16
wilcox.test(male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC5'], male_G46_fbias_rpkm$log2RPKM[male_G46_fbias_rpkm$cutoff=='FC7'])#W = 23165, p-value < 2.2e-16

wilcox.test(female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC1'], female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC3'])#W = 720900, p-value = 0.453
wilcox.test(female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC3'], female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC5'])#W = 27772, p-value = 0.005362
wilcox.test(female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC5'], female_G46_fbias_rpkm$log2RPKM[female_G46_fbias_rpkm$cutoff=='FC7'])#W = 10387, p-value = 0.0001314

###stats
g46_fbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='female')
g46_mbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='male')
g46_unbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='unbias')

wilcox.test(g46_fbias_rpkm$log2RPKM[g46_fbias_rpkm$Sex=='meanmale'], g46_fbias_rpkm$log2RPKM[g46_fbias_rpkm$Sex=='meanfemale'])#W = 3679200, p-value < 2.2e-16
wilcox.test(g46_mbias_rpkm$log2RPKM[g46_mbias_rpkm$Sex=='meanmale'], g46_mbias_rpkm$log2RPKM[g46_mbias_rpkm$Sex=='meanfemale'])#W = 1700900, p-value < 2.2e-16
wilcox.test(g46_unbias_rpkm$log2RPKM[g46_unbias_rpkm$Sex=='meanmale'], g46_unbias_rpkm$log2RPKM[g46_unbias_rpkm$Sex=='meanfemale'])#W = 184870000, p-value = 1.16e-09

### load data from gonad tissues
###sex bias cutoff plots
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amgonad/gonad_rpkm_FDR05_newcutoff_sorted.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XYtestis.XXovary<=-1 & kdata$FDR.XYtestis.XXovary< 0.05)
kdata2 <- subset(kdata, kdata$logFC.XYtestis.XXovary>=1 & kdata$FDR.XYtestis.XXovary< 0.05)

kdata1$meanXY <- (kdata1$A15MT1 + kdata1$A17MT1 + kdata1$A6MT1 + kdata1$A8MT1 + kdata1$A12MT1) /5 
kdata1$meanXX <- (kdata1$A10FO1 + kdata1$A17FO1 + kdata1$A2FO2 + kdata1$A6FO1 + kdata1$A8FO2) / 5 

kdata1_sub <- subset(kdata1, kdata1$meanXY>1)

kdata1 <- kdata1_sub

kdata1$log2meanXY <- log2(kdata1$meanXY)
kdata1$group1="meanmale"

kdata1_sub <- subset(kdata1, kdata1$meanXX>1)
kdata1 <- kdata1_sub

kdata1$log2meanXX <- log2(kdata1$meanXX)
kdata1$group2="meanfemale"

kdata2$meanXY <- (kdata2$A15MT1 + kdata2$A17MT1 + kdata2$A6MT1 + kdata2$A8MT1 + kdata2$A12MT1) /5 
kdata2$meanXX <- (kdata2$A10FO1 + kdata2$A17FO1 + kdata2$A2FO2 + kdata2$A6FO1 + kdata2$A8FO2) / 5 

kdata2_sub <- subset(kdata2,kdata2$meanXY>1 & kdata2$meanXX>1)
kdata2 <- kdata2_sub

kdata2$log2meanXY <- log2(kdata2$meanXY)
kdata2$group1="meanmale"
kdata2$log2meanXX <- log2(kdata2$meanXX)
kdata2$group2="meanfemale"

#write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias_cutoff_rpkm.txt", sep="\t", col.names=T)
#write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias_cutoff_rpkm.txt", sep="\t", col.names=T)

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias_cutoff_rpkm_rmslimit.txt", sep="\t", col.names=T)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias_cutoff_rpkm_rmslimit.txt", sep="\t", col.names=T)

###the version including sex-limited genes
gonad_mbias_rpkm_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias_cutoff_all_newfi.txt", header = TRUE)
str(gonad_mbias_rpkm_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/mbias_rpkm_gonad_cutoff.pdf")

p_gonadm <- ggplot(data = gonad_mbias_rpkm_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-15,-10,-5,0,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

gonad_fbias_rpkm_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias_rpkm_cutoff_newfi.txt", header = TRUE)
str(gonad_fbias_rpkm_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/fbias_rpkm_gonad_cutoff_newfi.pdf")

p_gonadf <- ggplot(data = gonad_fbias_rpkm_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-15,-10,-5,0,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

###stats
mbias_cutoff_ovary <- subset(gonad_mbias_rpkm_cutoff,gonad_mbias_rpkm_cutoff$sex=='meanfemale')
mbias_cutoff_testis <- subset(gonad_mbias_rpkm_cutoff,gonad_mbias_rpkm_cutoff$sex=='meanmale')

wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC1'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC3'])#W = 3155500, p-value < 2.2e-16
wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC3'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC5'])#W = 157720, p-value < 2.2e-16
wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC5'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC7'])#W = 16063, p-value = 8.454e-13

wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC1'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC3'])#W = 2167000, p-value = 0.03515
wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC3'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC5'])#W = 94164, p-value = 0.003794
wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC5'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC7'])#W = 8526, p-value = 0.005068

fbias_cutoff_ovary <- subset(gonad_fbias_rpkm_cutoff,gonad_fbias_rpkm_cutoff$sex=='meanfemale')
fbias_cutoff_testis <- subset(gonad_fbias_rpkm_cutoff,gonad_fbias_rpkm_cutoff$sex=='meanmale')

wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC1'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC3'])#W = 1600100, p-value = 3.807e-07
wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC3'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC5'])#W = 98967, p-value = 0.02975
wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC5'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC7']) #W = 30096, p-value = 5.585e-10

wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC1'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC3'])#W = 2251200, p-value < 2.2e-16
wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC3'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC5'])#W = 155170, p-value < 2.2e-16
wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC5'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC7'])#W = 57644, p-value = 3.37e-13

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/cutoff_rpkm_gonad.pdf")

ggarrange(p_gonadf, p_gonadm, labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

###the version removing sex-limited genes
gonad_mbias_rpkm_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias_cutoff_rpkm_rmslimit_fi.txt", header = TRUE)
str(gonad_mbias_rpkm_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/mbias_rpkm_gonad_cutoff_rmsexlimit.pdf")

p_gonadm <- ggplot(data = gonad_mbias_rpkm_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-0.5,0,2,5,10,15)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

gonad_fbias_rpkm_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias_cutoff_rpkm_rmslimit_fi.txt", header = TRUE)
str(gonad_fbias_rpkm_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/fbias_rpkm_gonad_cutoff_newfi_rmsexlimit.pdf")

p_gonadf <- ggplot(data = gonad_fbias_rpkm_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-0.5,0,2,5,10,15)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/cutoff_rpkm_gonad_rmsexlimit.pdf")

ggarrange(p_gonadf, p_gonadm, labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

###stats
mbias_cutoff_ovary <- subset(gonad_mbias_rpkm_cutoff,gonad_mbias_rpkm_cutoff$sex=='meanfemale')
mbias_cutoff_testis <- subset(gonad_mbias_rpkm_cutoff,gonad_mbias_rpkm_cutoff$sex=='meanmale')

wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC1'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC3'])#W = W = 539260, p-value = 3.498e-13
wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC3'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC5'])#W = 3384, p-value = 0.5552
wilcox.test(mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC5'], mbias_cutoff_ovary$log2RPKM[mbias_cutoff_ovary$cutoff=='FC7'])#W = 26, p-value = 0.2233

wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC1'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC3'])#***
wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC3'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC5'])#***
wilcox.test(mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC5'], mbias_cutoff_testis$log2RPKM[mbias_cutoff_testis$cutoff=='FC7'])#**

fbias_cutoff_ovary <- subset(gonad_fbias_rpkm_cutoff,gonad_fbias_rpkm_cutoff$sex=='meanfemale')
fbias_cutoff_testis <- subset(gonad_fbias_rpkm_cutoff,gonad_fbias_rpkm_cutoff$sex=='meanmale')

wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC1'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC3'])#***
wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC3'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC5'])#***
wilcox.test(fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC5'], fbias_cutoff_ovary$log2RPKM[fbias_cutoff_ovary$cutoff=='FC7']) #***

wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC1'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC3'])#***
wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC3'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC5'])#.
wilcox.test(fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC5'], fbias_cutoff_testis$log2RPKM[fbias_cutoff_testis$cutoff=='FC7'])#*

### load data from liver tissue
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amliver/rpkm_amliver_sorted_fi.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XYliver.XXliver<=-1)
kdata2 <- subset(kdata, kdata$logFC.XYliver.XXliver>=1)
kdata3 <- subset(kdata, kdata$logFC.XYliver.XXliver>-1 & kdata$logFC.XYliver.XXliver<1)

kdata1$meanXY <- (kdata1$logA12ML1 + kdata1$logA15ML1 + kdata1$logA17ML1 + kdata1$logA8ML1) /4 
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$logA12ML1 + kdata2$logA15ML1 + kdata2$logA17ML1 + kdata2$logA8ML1) /4 
kdata2$group1="meanmale"

kdata3$meanXY <- (kdata3$logA12ML1 + kdata3$logA15ML1 + kdata3$logA17ML1 + kdata3$logA8ML1) /4 
kdata3$group1="meanmale"

kdata1$meanXX <- (kdata1$logA12FL1 + kdata1$logA17FL1 + kdata1$logA2FL1 + kdata1$logA7FL1 + kdata1$logA8FL1) / 5 
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$logA12FL1 + kdata2$logA17FL1 + kdata2$logA2FL1 + kdata2$logA7FL1 + kdata2$logA8FL1) / 5 #XX female group at gonad
kdata2$group2="meanfemale"

kdata3$meanXX <- (kdata3$logA12FL1 + kdata3$logA17FL1 + kdata3$logA2FL1 + kdata3$logA7FL1 + kdata3$logA8FL1) / 5
kdata3$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_fbias.txt", sep="\t", col.names=F)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_mbias.txt", sep="\t", col.names=F)
write.table(kdata3, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_unbias.txt", sep="\t", col.names=F)

liver_fbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_fbias1.txt", header=TRUE)
liver_fbias_all$bias = "female"
liver_mbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_mbias1.txt", header=TRUE)
liver_mbias_all$bias = "male"
liver_unbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_unbias1.txt", header=TRUE)
liver_unbias_all$bias = "unbias"

write.table(liver_fbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_fbias_all.txt", sep="\t", col.names=F)
write.table(liver_mbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_mbias_all.txt", sep="\t", col.names=F)
write.table(liver_unbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_unbias_all.txt", sep="\t", col.names=F)

liver_sbunbias_rpkm <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_sbunbias_all_fi.txt", header = TRUE)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbunbias_logrpkm_liver.pdf")

sts <- boxplot.stats(liver_sbunbias_rpkm$log2RPKM)$stats

p4 <- ggplot(data = liver_sbunbias_rpkm, aes(x=Bias, y= log2RPKM, fill=Sex)) + 
  geom_boxplot(notch=TRUE,position=position_dodge(0.8))  +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  coord_cartesian(ylim = c(sts[1]*1.5,max(sts)*1.5)) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
#  theme(legend.position="none")
  
dev.off()

####
liver_fbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='female')
liver_mbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='male')
liver_unbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='unbias')

wilcox.test(liver_fbias_rpkm$log2RPKM[liver_fbias_rpkm$Sex=='meanmale'], liver_fbias_rpkm$log2RPKM[liver_fbias_rpkm$Sex=='meanfemale'])#W = 282800, p-value < 2.2e-16
wilcox.test(liver_mbias_rpkm$log2RPKM[liver_mbias_rpkm$Sex=='meanmale'], liver_mbias_rpkm$log2RPKM[liver_mbias_rpkm$Sex=='meanfemale']) #W = 285950, p-value < 2.2e-16
wilcox.test(liver_unbias_rpkm$log2RPKM[liver_unbias_rpkm$Sex=='meanmale'],liver_unbias_rpkm$log2RPKM[liver_unbias_rpkm$Sex=='meanfemale']) #W = 141100000, p-value = 0.319

###sex-bias cutoff figures

kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amliver/rpkm_amliver_sorted_cutoff_fi.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XYliver.XXliver<=-1 & kdata$FDR.XYliver.XXliver < 0.05)
kdata2 <- subset(kdata, kdata$logFC.XYliver.XXliver>=1 & kdata$FDR.XYliver.XXliver < 0.05)

kdata1$meanXY <- (kdata1$A12ML1 + kdata1$A15ML1 + kdata1$A17ML1 + kdata1$A8ML1) /4 
kdata1$log2meanXY <- log2(kdata1$meanXY + 0.0001)
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$A12ML1 + kdata2$A15ML1 + kdata2$A17ML1 + kdata2$A8ML1) /4 
kdata2$log2meanXY <- log2(kdata2$meanXY + 0.0001)
kdata2$group1="meanmale"


kdata1$meanXX <- (kdata1$A12FL1 + kdata1$A17FL1 + kdata1$A2FL1 + kdata1$A7FL1 + kdata1$A8FL1) / 5 
kdata1$log2meanXX <- log2(kdata1$meanXX + 0.0001)
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$A12FL1 + kdata2$A17FL1 + kdata2$A2FL1 + kdata2$A7FL1 + kdata2$A8FL1) / 5 #XX female group at gonad
kdata2$log2meanXX <- log2(kdata2$meanXX + 0.0001)
kdata2$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_fbias2_cutoff2.txt", sep="\t", col.names=T)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_mbias2_cutoff2.txt", sep="\t", col.names=T)

liver_mbias_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_mbias2_cutoff_inclusexlimit_fi.txt", header=TRUE)
str(liver_mbias_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_liver_mbias.pdf")

p1 <-  ggplot(data =liver_mbias_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-20,-10, -5, 0,2,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()


liver_fbias_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/liver_sb_cutoff_inclusexlimit_fi.txt", header=TRUE)
str(liver_fbias_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_liver_fbias.pdf")

p2 <- ggplot(data =liver_fbias_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-20,-10, -5, 0,2,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_liver_inclusexlimit.pdf")

ggarrange(p2, p1, labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

###stats
liver_fbias_cutoff_male <- subset(liver_fbias_cutoff, liver_fbias_cutoff$sex=='meanmale')
liver_fbias_cutoff_female  <- subset(liver_fbias_cutoff , liver_fbias_cutoff $sex=='meanfemale')

wilcox.test(liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC1'], liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC3'])#W = 1997, p-value = 3.833e-05
wilcox.test(liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC3'], liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC5'])#W = 263, p-value = 0.2815
wilcox.test(liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC5'], liver_fbias_cutoff_male$log2RPKM[liver_fbias_cutoff_male$cutoff=='FC7'])#W = 171, p-value = 0.01125

wilcox.test(liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC1'], liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC3'])#W = 1544, p-value = 0.2162
wilcox.test(liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC3'], liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC5'])#W = 130, p-value = 0.0002754
wilcox.test(liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC5'], liver_fbias_cutoff_female$log2RPKM[liver_fbias_cutoff_female$cutoff=='FC7'])#W = 127, p-value = 0.5088

liver_mbias_cutoff_male <- subset(liver_mbias_cutoff, liver_mbias_cutoff$sex=='meanmale')
liver_mbias_cutoff_female  <- subset(liver_mbias_cutoff , liver_mbias_cutoff $sex=='meanfemale')

wilcox.test(liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC1'], liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC3'])#W = 705, p-value = 0.003533
wilcox.test(liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC3'], liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC5'])#W = 13, p-value = 0.2799
wilcox.test(liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC5'], liver_mbias_cutoff_male$log2RPKM[liver_mbias_cutoff_male$cutoff=='FC7'])#W = 6, p-value = 0.5333

wilcox.test(liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC1'], liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC3'])#W = 828, p-value = 1.278e-05
wilcox.test(liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC3'], liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC5'])#W = 25, p-value = 0.7531
wilcox.test(liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC5'], liver_mbias_cutoff_female$log2RPKM[liver_mbias_cutoff_female$cutoff=='FC7'])#W = 8, p-value = 0.1333

###load G43 datasets
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg43/rpkm_ammg43_fi.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
kdata2 <- subset(kdata, kdata$logFC.XY43.XX43>=1)
kdata3 <- subset(kdata, kdata$logFC.XY43.XX43>-1 & kdata$logFC.XY43.XX43<1)

kdata1$meanXY <- (kdata1$logAm1_434 + kdata1$logAm2_434 + kdata1$logAm4_434) /3 
 kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$logAm1_434 + kdata2$logAm2_434 + kdata2$logAm4_434) /3 
kdata2$group1="meanmale"

kdata3$meanXY <- (kdata3$logAm1_434 + kdata3$logAm2_434 + kdata3$logAm4_434) /3 
kdata3$group1="meanmale"

kdata1$meanXX <- (kdata1$logAm2_433 + kdata1$logAm4_435 + kdata1$logAm5_433)/3
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$logAm2_433 + kdata2$logAm4_435 + kdata2$logAm5_433)/3  #XX female group at gonad
kdata2$group2="meanfemale"

kdata3$meanXX <- (kdata3$logAm2_433 + kdata3$logAm4_435 + kdata3$logAm5_433)/3
kdata3$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_fbias.txt", sep="\t", col.names=F)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_mbias.txt", sep="\t", col.names=F)
write.table(kdata3, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_unbias.txt", sep="\t", col.names=F)

g43_fbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_fbias1.txt", header=TRUE)
g43_fbias_all$bias = "female"
g43_mbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_mbias1.txt", header=TRUE)
g43_mbias_all$bias = "male"
g43_unbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_unbias1.txt", header=TRUE)
g43_unbias_all$bias = "unbias"

write.table(g43_fbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_fbias_all.txt", sep="\t", col.names=F)
write.table(g43_mbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_mbias_all.txt", sep="\t", col.names=F)
write.table(g43_unbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_unbias_all.txt", sep="\t", col.names=F)

g43_sbunbias_rpkm <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_sbunbias_rpkm_fi.txt", header = TRUE)
str(g43_sbunbias_rpkm)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbunbias_logrpkm_g43.pdf")

sts <- boxplot.stats(g43_sbunbias_rpkm$log2RPKM)$stats

p1 <- ggplot(data = g43_sbunbias_rpkm, aes(x=Bias, y= log2RPKM, fill=Sex)) + 
  geom_boxplot(notch=TRUE,outlier.colour = NA, position=position_dodge(0.8))  +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  coord_cartesian(ylim = c(sts[1]*1.8,max(sts)*1.3)) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

####
g43_fbias_rpkm <- subset(g43_sbunbias_rpkm, g43_sbunbias_rpkm$Bias=='female')
g43_mbias_rpkm <- subset(g43_sbunbias_rpkm, g43_sbunbias_rpkm$Bias=='male')
g43_unbias_rpkm <- subset(g43_sbunbias_rpkm, g43_sbunbias_rpkm$Bias=='unbias')

wilcox.test(g43_fbias_rpkm$log2RPKM[g43_fbias_rpkm$Sex=='meanmale'], g43_fbias_rpkm$log2RPKM[g43_fbias_rpkm$Sex=='meanfemale'])#W = 99431, p-value < 2.2e-16
wilcox.test(g43_mbias_rpkm$log2RPKM[g43_mbias_rpkm$Sex=='meanmale'], g43_mbias_rpkm$log2RPKM[g43_mbias_rpkm$Sex=='meanfemale']) #W = 245000, p-value < 2.2e-16
wilcox.test(g43_unbias_rpkm$log2RPKM[g43_unbias_rpkm$Sex=='meanmale'],g43_unbias_rpkm$log2RPKM[g43_unbias_rpkm$Sex=='meanfemale']) #W = 288400000, p-value = 0.3782

##put four figures in one page
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbunbias_logrpkm_g43G46gonadliver.pdf")

ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

dev.off()

###G43 for sex-bias cut off figures

kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg43/rpkm_ammg43_forcutoff.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XY43.XX43<=-1 & kdata$FDR.XY43.XX43 < 0.05)
kdata2 <- subset(kdata, kdata$logFC.XY43.XX43>=1 & kdata$FDR.XY43.XX43 < 0.05)

kdata1$meanXY <- (kdata1$Am1_434 + kdata1$Am2_434 + kdata1$Am4_434) /3 
kdata1$log2meanXY <- log2(kdata1$meanXY + 0.0001)
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$Am1_434 + kdata2$Am2_434 + kdata2$Am4_434) /3 
kdata2$log2meanXY <- log2(kdata2$meanXY + 0.0001)
kdata2$group1="meanmale"

kdata1$meanXX <- (kdata1$Am2_433 + kdata1$Am4_435 + kdata1$Am5_433)/3
kdata1$log2meanXX <- log2(kdata1$meanXX + 0.0001)
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$Am2_433 + kdata2$Am4_435 + kdata2$Am5_433)/3  #XX female group at gonad
kdata2$log2meanXX <- log2(kdata2$meanXX + 0.0001)
kdata2$group2="meanfemale"


write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_fbias2_cutoff.txt", sep="\t", col.names=T)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_mbias2_cutoff.txt", sep="\t", col.names=T)

g43_fbias_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_fbias2_cutoff_inclsexlimit_fi.txt", header=TRUE)
str(g43_fbias_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_g43_fbias.pdf")

pg43_f <- ggplot(data = g43_fbias_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-15,-10, -5, 0,2,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

g43_mbias_cutoff <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g43_mbias_cutoff_inclusexlimit_fi.txt", header=TRUE)
str(g43_mbias_cutoff)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_g43_mbias.pdf")

pg43_m <- ggplot(data = g43_mbias_cutoff, aes(x=cutoff, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0), width=0.9,alpha=0.9) +
  coord_cartesian(ylim = c(-15,-10, -5, 0,2,5,10)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

dev.off()

####
g43_fbias_cutoff_male <- subset(g43_mbias_cutoff, g43_mbias_cutoff$sex=='meanmale')
g43_mbias_cutoff_female <- subset(g43_mbias_cutoff, g43_mbias_cutoff$sex=='meanfemale')

wilcox.test(g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC1'], g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC3'])#w== 57, p-value = 0.3562
wilcox.test(g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC3'], g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC5'])#W = 14, p-value = 0.9371
wilcox.test(g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC5'], g43_mbias_cutoff_male$log2RPKM[g43_mbias_cutoff_male$cutoff=='FC7'])#W = 15, p-value = 0.6303

wilcox.test(g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC1'], g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC3'])#W = 66, p-value = 0.09472
wilcox.test(g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC3'], g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC5'])#W = 24, p-value = 0.1608
wilcox.test(g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC5'], g43_mbias_cutoff_female$log2RPKM[g43_mbias_cutoff_female$cutoff=='FC7'])#W = 23, p-value = 0.01943

g43_fbias_cutoff_male <- subset(g43_fbias_cutoff, g43_fbias_cutoff$sex=='meanmale')
g43_fbias_cutoff_female <- subset(g43_fbias_cutoff, g43_fbias_cutoff$sex=='meanfemale')

wilcox.test(g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC1'], g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC3'])#W = 788, p-value = 0.5908
wilcox.test(g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC3'], g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC5'])#W = 1416, p-value = 0.002866
wilcox.test(g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC5'], g43_fbias_cutoff_female$log2RPKM[g43_fbias_cutoff_female$cutoff=='FC7'])#W = 307, p-value = 0.01566

wilcox.test(g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC1'], g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC3'])#W = 1132, p-value = 0.0001738
wilcox.test(g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC3'], g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC5'])#W = 2656, p-value = 0.003664
wilcox.test(g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC5'], g43_fbias_cutoff_male$log2RPKM[g43_fbias_cutoff_male$cutoff=='FC7'])#W = 392, p-value = 1.881e-05

##put four figures in one page
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbcutoff_logrpkm_g43_inclusexlimit.pdf")

ggarrange(pg43_f,pg43_m, labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()
