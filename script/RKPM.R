# laading R libraries
library(edgeR)
library(gplots)
library(ggplot2)
install.packages("gdata", dependencies=TRUE) #for cbind of multiple data frames
library(gdata)
install.packages("plyr")
library(plyr)

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
p + geom_violin(scale="width",fill="grey", trim=F) + 
  geom_hline(yintercept=-1,linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") +
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


col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/rpkm_list_count_sorted.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
kdata1 <- subset(kdata, kdata$logFC.XY46.XX46<=-1)
kdata2 <- subset(kdata, kdata$logFC.XY46.XX46>=1)
kdata3 <- subset(kdata, kdata$logFC.XY46.XX46>-1 & kdata$logFC.XY46.XX46<1)

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
ggplot(data = g46_bias_all, aes(x=Bias, y= log2RPKM, fill=Sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA) +
  coord_cartesian(ylim = c(sts[1]*1.3,max(sts)*1.3)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12))
dev.off()

###stats
g46_fbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='female')
g46_mbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='male')
g46_unbias_rpkm <- subset(g46_bias_all, g46_bias_all$Bias=='unbias')

wilcox.test(g46_fbias_rpkm$log2RPKM[g46_fbias_rpkm$Sex=='meanmale'], g46_fbias_rpkm$log2RPKM[g46_fbias_rpkm$Sex=='meanfemale'])#W = 3679200, p-value < 2.2e-16
wilcox.test(g46_mbias_rpkm$log2RPKM[g46_mbias_rpkm$Sex=='meanmale'], g46_mbias_rpkm$log2RPKM[g46_mbias_rpkm$Sex=='meanfemale'])#W = 1700900, p-value < 2.2e-16
wilcox.test(g46_unbias_rpkm$log2RPKM[g46_unbias_rpkm$Sex=='meanmale'], g46_unbias_rpkm$log2RPKM[g46_unbias_rpkm$Sex=='meanfemale'])#W = 184870000, p-value = 1.16e-09

### load data from gonad tissues
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amgonad/rpkm_amm_gonad.txt", header=TRUE)
str(kdata)

kdata1 <- subset(kdata, kdata$logFC.XYtestis.XXovary<=-1)
kdata2 <- subset(kdata, kdata$logFC.XYtestis.XXovary>=1)
kdata3 <- subset(kdata, kdata$logFC.XYtestis.XXovary>-1 & kdata$logFC.XYtestis.XXovary<1)

kdata1$meanXY <- (kdata1$logA15MT1 + kdata1$logA17MT1 + kdata1$logA6MT1 + kdata1$logA8MT1 + kdata1$logA12MT1) /5 
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$logA15MT1 + kdata2$logA17MT1 + kdata2$logA6MT1 + kdata2$logA8MT1 + kdata2$logA12MT1) / 5 
kdata2$group1="meanmale"

kdata3$meanXY <- (kdata3$logA15MT1 + kdata3$logA17MT1 + kdata3$logA6MT1 + kdata3$logA8MT1 + kdata3$logA12MT1) /5
kdata3$group1="meanmale"

kdata1$meanXX <- (kdata1$logA10FO1 + kdata1$logA17FO1 + kdata1$logA2FO2 + kdata1$logA6FO1 + kdata1$logA8FO2) / 5 
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$logA10FO1 + kdata2$logA17FO1 + kdata2$logA2FO2 + kdata2$logA6FO1 + kdata2$logA8FO2) / 5 #XX female group at gonad
kdata2$group2="meanfemale"

kdata3$meanXX <- (kdata3$logA10FO1 + kdata3$logA17FO1 + kdata3$logA2FO2 + kdata3$logA6FO1 + kdata3$logA8FO2) / 5
kdata3$group2="meanfemale"

write.table(kdata1, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias.txt", sep="\t", col.names=F)
write.table(kdata2, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias.txt", sep="\t", col.names=F)
write.table(kdata3, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_unbias.txt", sep="\t", col.names=F)

gonad_fbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias2.txt", header=TRUE)
gonad_fbias_all$bias = "female"
gonad_mbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias2.txt", header=TRUE)
gonad_mbias_all$bias = "male"
gonad_unbias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_unbias2.txt", header=TRUE)
gonad_unbias_all$bias = "unbias"

write.table(gonad_fbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_fbias_all.txt", sep="\t", col.names=F)
write.table(gonad_mbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_mbias_all.txt", sep="\t", col.names=F)
write.table(gonad_unbias_all, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_unbias_all.txt", sep="\t", col.names=F)

gonad_all_rpkm <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/gonad_sbunbias_rpkm_fi.txt", header = TRUE)
str(gonad_all_rpkm)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sbunbias_logrpkm_gonad.pdf")

sts <- boxplot.stats(gonad_all_rpkm$log2RPKM)$stats
ggplot(data = gonad_all_rpkm, aes(x=bias, y= log2RPKM, fill=sex)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA) +
  coord_cartesian(ylim = c(sts[1]*1.3,max(sts)*1.3)) +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  theme_set(theme_bw(base_size=12))
dev.off()

###stats
gonad_fbias_rpkm <- subset(gonad_all_rpkm, gonad_all_rpkm$bias=='female')
gonad_mbias_rpkm <- subset(gonad_all_rpkm, gonad_all_rpkm$bias=='male')
gonad_unbias_rpkm <- subset(gonad_all_rpkm, gonad_all_rpkm$bias=='unbias')

wilcox.test(gonad_fbias_rpkm$log2RPKM[gonad_fbias_rpkm$sex=='meanmale'], gonad_fbias_rpkm$log2RPKM[gonad_fbias_rpkm$sex=='meanfemale'])#W = 6990300, p-value < 2.2e-16
wilcox.test(gonad_mbias_rpkm$log2RPKM[gonad_mbias_rpkm$sex=='meanmale'], gonad_mbias_rpkm$log2RPKM[gonad_mbias_rpkm$sex=='meanfemale']) #W = 24259000, p-value < 2.2e-16
wilcox.test(gonad_unbias_rpkm$log2RPKM[gonad_unbias_rpkm$sex=='meanmale'],gonad_unbias_rpkm$log2RPKM[gonad_unbias_rpkm$sex=='meanfemale']) #W = 56847000, p-value = 0.4779

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

ggplot(data = liver_sbunbias_rpkm, aes(x=Bias, y= log2RPKM, fill=Sex)) + 
  geom_boxplot(notch=TRUE,outlier.colour = NA)  +
  scale_fill_manual(values = c("red3","dodgerblue3")) +
  coord_cartesian(ylim = c(sts[1]*1.5,max(sts)*1.5)) +
  theme_set(theme_bw(base_size=12))
  
dev.off()

####
liver_fbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='female')
liver_mbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='male')
liver_unbias_rpkm <- subset(liver_sbunbias_rpkm, liver_sbunbias_rpkm$Bias=='unbias')

wilcox.test(liver_fbias_rpkm$log2RPKM[liver_fbias_rpkm$Sex=='meanmale'], liver_fbias_rpkm$log2RPKM[liver_fbias_rpkm$Sex=='meanfemale'])#W = 282800, p-value < 2.2e-16
wilcox.test(liver_mbias_rpkm$log2RPKM[liver_mbias_rpkm$Sex=='meanmale'], liver_mbias_rpkm$log2RPKM[liver_mbias_rpkm$Sex=='meanfemale']) #W = 285950, p-value < 2.2e-16
wilcox.test(liver_unbias_rpkm$log2RPKM[liver_unbias_rpkm$Sex=='meanmale'],liver_unbias_rpkm$log2RPKM[liver_unbias_rpkm$Sex=='meanfemale']) #W = 141100000, p-value = 0.319


###load G43 datasets
