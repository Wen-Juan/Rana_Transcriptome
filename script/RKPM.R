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
p + geom_violin(scale="width",fill="white", trim=F) + 
  geom_hline(yintercept=-1,linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") +
  #+ geom_hline(yintercept=-2,linetype="dashed") + geom_hline(yintercept=2, linetype="dashed") + 
  theme(panel.background = element_rect(fill='blue')) + scale_y_continuous(breaks=c(-15,-10, -5, -1, 0, 1, 5, 10, 15)) +theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
dev.off()

### dataset from G46


col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg46nosr/LogCPM_0.05_Amg46nosr copy.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
kdata1 <- subset(kdata, kdata$logFC.XY46.XX46<=-1)
kdata2 <- subset(kdata, kdata$logFC.XY46.XX46>=1 & kdata$logFC.XY46.XX46>=1)
kdata3 <- subset(kdata, kdata$logFC.XY46.XX46>-1 & kdata$logFC.XY46.XX46<1)

kdata1$meanXY <- (kdata1$Am2_463 + kdata1$Am4_461 + kdata1$Am6_462) / 3 #XY male group at G46
kdata1$group1="meanmale"

kdata2$meanXY <- (kdata2$Am2_463 + kdata2$Am4_461 + kdata2$Am6_462) / 3 #XY male group at G46
kdata2$group1="meanmale"

kdata3$meanXY <- (kdata3$Am2_463 + kdata3$Am4_461 + kdata3$Am6_462) / 3 #XY male group at G46
kdata3$group1="meanmale"
#group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY male group at G46
#group1.expr.minus <- (map.data$Am1_434 + map.data$Am2_434 + map.data$Am4_434) / 3 #XY male group at G43
#group1.expr.minus <- (map.data$A15MT1 + map.data$A17MT1 + map.data$A6MT1 + map.data$A8MT1 + map.data$A12MT1) / 5 #XY male group in gonad
#group1.expr.minus <- (map.data$A10MB + map.data$A15MB + map.data$A16MB + map.data$A17MB + map.data$A8MB) /5 #XY in brain
#group1.expr.minus <- (map.data$A12ML1 + map.data$A17ML1 + map.data$A15ML1 +map.data$A8ML1)/4 #XY linver

kdata1$meanXX <- (kdata1$Am2_464+kdata1$Am6_464)/2 #XX male group at G46
kdata1$group2="meanfemale"

kdata2$meanXX <- (kdata2$Am2_464+kdata2$Am6_464)/2 #XX male group at G46
kdata2$group2="meanfemale"

kdata3$meanXX <- (kdata3$Am2_464+kdata3$Am6_464)/2 #XY male group at G46
kdata3$group2="meanfemale"
#group2.expr.minus <- (map.data$Am2_433 + map.data$Am4_435 + map.data$Am5_433) / 3 #XX female group at G43
#group2.expr.minus <- (map.data$A10FO1 + map.data$A17FO1 + map.data$A2FO2 + map.data$A8FO2 + map.data$A6FO1) / 5 #XY male group in gonad
#group2.expr.minus <- (map.data$A10FB + map.data$A12FB + map.data$A15FB + map.data$A16FB + map.data$A17FB) /5 #XX in brain
#group2.expr.minus <- (map.data$A12FL1 + map.data$A17FL1 + map.data$A2FL1 + map.data$A7FL1 +map.data$A8FL1)/5 #XX linver

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

g46_bias_all <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/CPM/g46_all_bias.txt", header=TRUE)

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sexbias_logcpmg46_alltog.pdf")

ggplot(data = g46_bias_all, aes(x=bias, y= LogCPM, fill=sex)) + 
  geom_boxplot(notch=TRUE)  +
  scale_fill_manual(values = c("red","blue")) +
  theme(axis.title.x = element_text(size=15,colour = "black"),axis.title.y = element_text(size=15,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sexbias_logcpmg46_fbias.pdf")

ggplot(g46_fbias_all, aes(x=sex, y=logCPM, fill=sex)) + 
  scale_fill_manual(values = c("red","blue")) +
  theme(legend.position="none") +
  geom_boxplot(notch=TRUE) +
  ylim(-5,13) +
  labs(y='LogCPM') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sexbias_logcpmg46_mbias.pdf")

ggplot(g46_mbias_all, aes(x=sex, y=logCPM, fill=sex)) + 
  scale_fill_manual(values = c("red","blue")) +
  theme(legend.position="none") +
  geom_boxplot(notch=TRUE) +
  ylim(-5,13) +
  labs(y='LogCPM') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/sexbias_logcpmg46_unbias.pdf")

ggplot(g46_unbias_all, aes(x=sex, y=logCPM, fill=sex)) + 
  scale_fill_manual(values = c("red","blue")) +
  theme(legend.position="none") +
  geom_boxplot(notch=TRUE) +
  ylim(-5,13) +
  labs(y='LogCPM') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()
