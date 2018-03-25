#intall R packages and load libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
library("ggplot2")


#set working directory
datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/'
results <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/'

#loading data
Amm_all_TPM <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/Amm_all_TPM.txt",header = T)
str(Amm_all_TPM)

###tissue number
######developmental stages, 7 tissues in total: G23, G27, G31, G43_female(G43_XX), G43_male(G43_XY), G46_female(G46_XX), G46_male(G46_XY) 
######adult tissues, 4 tissues in total: testis, ovary, brain, liver
#for the new tau calculation, we split each tissue category in female and male tissues, so in total 16 tissues.

###function requires data frame with expression values
#mean values between replicates are calculated.

Amm_all_TPM$G23_fmean <- (Amm_all_TPM$G23_F.1 +Amm_all_TPM$G23_F.2+Amm_all_TPM$G23_F)/3
Amm_all_TPM$G23_mmean <- (Amm_all_TPM$G23_M.1 +Amm_all_TPM$G23_M.2+Amm_all_TPM$G23_M)/3

Amm_all_TPM$G27_fmean <- (Amm_all_TPM$G27_F.1 +Amm_all_TPM$G27_F.2+Amm_all_TPM$G27_F)/3
Amm_all_TPM$G27_mmean <- (Amm_all_TPM$G27_M.1 +Amm_all_TPM$G27_M.2+Amm_all_TPM$G27_M)/3

Amm_all_TPM$G31_fmean <- (Amm_all_TPM$G31_F.1 +Amm_all_TPM$G31_F.2+Amm_all_TPM$G31_F.3+Amm_all_TPM$G31_F)/4
Amm_all_TPM$G31_mmean <- (Amm_all_TPM$G31_M.1 +Amm_all_TPM$G31_M.2+Amm_all_TPM$G31_M.3+Amm_all_TPM$G31_M.4+Amm_all_TPM$G31_M)/5

Amm_all_TPM$G43_M_mean <- (Amm_all_TPM$G43_M + Amm_all_TPM$G43_M.1 + Amm_all_TPM$G43_M.2)/3
Amm_all_TPM$G43_F_mean <- (Amm_all_TPM$G43_F + Amm_all_TPM$G43_F.1 + Amm_all_TPM$G43_F.2)/3
Amm_all_TPM$G46_F_mean <- (Amm_all_TPM$G46_F + Amm_all_TPM$G46_F.1)/2
Amm_all_TPM$G46_M_mean <- (Amm_all_TPM$G46_M + Amm_all_TPM$G46_M.1 + Amm_all_TPM$G46_M.2)/3

Amm_all_TPM$Testis_mean <- (Amm_all_TPM$Testis + Amm_all_TPM$Testis.1 + Amm_all_TPM$Testis.2 +Amm_all_TPM$Testis.3 +Amm_all_TPM$Testis.4)/5
Amm_all_TPM$Ovary_mean <- (Amm_all_TPM$Ovary + Amm_all_TPM$Ovary.1 + Amm_all_TPM$Ovary.2 +Amm_all_TPM$Ovary.3 +Amm_all_TPM$Ovary.4)/5

Amm_all_TPM$Brain_fmean <- (Amm_all_TPM$Brain_F + Amm_all_TPM$Brain_F.1 + Amm_all_TPM$Brain_F.2 +Amm_all_TPM$Brain_F.3 +Amm_all_TPM$Brain_F.4 +Amm_all_TPM$Brain_F)/5
Amm_all_TPM$Brain_mmean <- (Amm_all_TPM$Brain_M + Amm_all_TPM$Brain_M.1 + Amm_all_TPM$Brain_M.2 +Amm_all_TPM$Brain_M.3 +Amm_all_TPM$Brain_M.4 +Amm_all_TPM$Brain_M)/5

Amm_all_TPM$Liver_fmean <- (Amm_all_TPM$Liver_F + Amm_all_TPM$Liver_F.1 + Amm_all_TPM$Liver_F.2 +Amm_all_TPM$Liver_F.3 +Amm_all_TPM$Liver_F.4 +Amm_all_TPM$Liver_F)/5
Amm_all_TPM$Liver_mmean <- (Amm_all_TPM$Liver_M + Amm_all_TPM$Liver_M.1 + Amm_all_TPM$Liver_M.2 +Amm_all_TPM$Liver_M.3 +Amm_all_TPM$Liver_M)/5

Amm_all_TPM$all_mean <- (Amm_all_TPM$G23_fmean+Amm_all_TPM$G23_mmean+Amm_all_TPM$G27_fmean+Amm_all_TPM$G27_mmean+Amm_all_TPM$G31_fmean+Amm_all_TPM$G31_mmean+Amm_all_TPM$G43_M_mean+Amm_all_TPM$G43_F_mean+Amm_all_TPM$G46_F_mean+Amm_all_TPM$G46_M_mean+Amm_all_TPM$Testis_mean+Amm_all_TPM$Ovary_mean+Amm_all_TPM$Liver_fmean+Amm_all_TPM$Liver_mmean+Amm_all_TPM$Brain_fmean+Amm_all_TPM$Brain_mmean)/16
str(Amm_all_TPM$all_mean)

#select TPM>=1 in at least one tissue out of 16 tissues, average TPM>=0.1
Amm_all_TPM <- Amm_all_TPM[Amm_all_TPM$all_mean>=0.1,]
a1 <- Amm_all_TPM
m <- 15
rownames(a1) <- a1[,1]
a1 <-a1[,-1]
a1 <- round(a1,3)
dim1 <- dim(a1)[1]
a2 <- a1[,c(62,63:78)]
dim1 <- dim(a2)[1]

s <-matrix(0, nrow=dim(a2)[1], ncol=2, dimnames = list(rownames(a2),c("Tissue_specific_index", "Tissue_name")))

#with log transformation, reduce effects of extreme data and variation inbetween.
for(i in 1:dim1){ 
  s[i,1] <- (sum(1 - log(1+(a2[i,]))/max(log(1+(a2[i,])))))/ m   
  s[i,2] <- names(which.max(log(1+(a2[i,]))))
}

write.table(s, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm.txt",row.names=TRUE, sep="\t", quote= FALSE)

#plot tau
tau_amm <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm.txt", header=T)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau//tau_amm_logtransf.pdf", width=8, height=8)
ggplot(tau_amm, aes(Tissue_specific_index)) +
  geom_histogram(binwidth = 0.01,fill=I("blue"),col=I("blue"))
dev.off()

#without log tansformation, data might be with higher variation
for(i in 1:dim1){ 
  s[i,1] <- (sum(1 - (a2[i,])/max(a2[i,])))/ m   
  s[i,2] <- names(which.max(a2[i,]))
}

write.table(s, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm_withoutlogtransf.txt",row.names=TRUE, sep="\t", quote= FALSE)

#make ggplot for non-log transformed tau distribution
tau_amm_nolog <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm_withoutlogtransf.txt", header=T)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm_withoutlogtransf.pdf", width=8, height=8)
ggplot(tau_amm_nolog, aes(Tissue_specific_index)) +
  geom_histogram(binwidth = 0.01,fill=I("blue"),col=I("blue"))
dev.off()

#subset the tissue specific genes (high tau values, >=0.9), and also genes with more pleotropic effects (tau<0.2).
s2 <- subset(tau_amm,tau_amm[,2] >= 0.9)
write.table(s2, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm_log_TSI0.9.txt",row.names=TRUE, sep="\t", quote= FALSE)
s3 <- subset(tau_amm,tau_amm[,2]<=0.2)
write.table(s3, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/tau_amm_log_TSI0.2.txt",row.names=TRUE, sep="\t", quote= FALSE)

##############
#load data from folder G46
##############

g46_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/G46/Amg46nosr_log1_sb_un_tau_fi.txt", header = TRUE)
str(g46_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_amm_g46.pdf", width=8, height=8)
ggplot(g46_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1.25) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()


y1 <- lm(sqrt(tau)~sqrt(abs(g46_tau$logFC.XY46.XX46))*bias, data=g46_tau)
anova(y1) 

#############
#############with unbias genes
Df Sum Sq Mean Sq  F value    Pr(>F)    
sqrt(abs(g46_tau$logFC.XY46.XX46))          1  65.92  65.923 3775.405 < 2.2e-16 ***
  bias                                        2   2.67   1.336   76.499 < 2.2e-16 ***
  sqrt(abs(g46_tau$logFC.XY46.XX46)):bias     2   2.23   1.114   63.779 < 2.2e-16 ***
  Residuals                               24400 426.05   0.017   
####

##investigate whether dnds is influences by tau or sex bias, or the interaction between the two.
g46_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/G46/Amg46nosr_log1_sb_un_tau_dnds_fi.txt", header = TRUE)
str(g46_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_amm_g46_dnds.pdf", width=8, height=8)
ggplot(g46_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='male'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 441010, p-value = 8.044e-08
wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='female'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 4095500, p-value = 0.1994
wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='male'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='female'],exact = FALSE) 
#W = 152450, p-value = 3.095e-06

y1 <- lm(sqrt(dNdS)~sqrt(tau)*bias, g46_tau_dnds)
anova(y1)

###########
Df Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)         1  1.941 1.94074 212.689 < 2.2e-16 ***
  bias              2  0.200 0.10001  10.960 1.770e-05 ***
  sqrt(tau):bias    2  0.375 0.18746  20.544 1.275e-09 ***
  Residuals      6589 60.123 0.00912    
###########

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g46_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC.XY46.XX46)", y="Tau", color ="bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g46_taudnds_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=g46_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g46_sb.pdf", width=8, height=8)
ggplot2.scatterplot(data=g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="LogFC.XY46.XX46)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g46_sbtau.pdf", width=8, height=8)
ggplot2.scatterplot(data=g46_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="dNdS)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

#stats
wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='unbias'],exact = FALSE) 
#W = 9537200, p-value < 2.2e-16

wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 1494500, p-value = 0.01684

wilcox.test(g46_tau$tau[g46_tau$bias=='unbias'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 27165000, p-value < 2.2e-16

##############
#load data from folder G43
##############
g43_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/G43/Amg43_sb_un_tau_sorted_fi.txt", header = TRUE)
str(g43_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/new_tau_amm_g43.pdf", width=8, height=8)
ggplot(g43_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

#wilcox test for tau
wilcox.test(g43_tau$tau[g43_tau$bias=='male'],g43_tau$tau[g43_tau$bias=='female'],exact = FALSE) 
#W = 456, p-value = 3.07e-12
wilcox.test(g43_tau$tau[g43_tau$bias=='unbias'],g43_tau$tau[g43_tau$bias=='female'],exact = FALSE) 
#W = 150230, p-value < 2.2e-16
wilcox.test(g43_tau$tau[g43_tau$bias=='unbias'],g43_tau$tau[g43_tau$bias=='male'],exact = FALSE) 
#W = 166880, p-value = 1.002e-07

####some stats

str(g43_tau)
qqnorm(sqrt(abs(g43_tau$logFC.XY43.XX43)))
hist(sqrt(abs(g43_tau$logFC.XY43.XX43)))
qqnorm(sqrt(abs(g43_tau$tau)))
hist(sqrt(abs(g43_tau$tau)))

y <- lm(sqrt(tau)~sqrt(abs(logFC.XY43.XX43))*bias, data=g43_tau)
anova(y)

######
Response: sqrt(tau)
Df Sum Sq Mean Sq   F value    Pr(>F)    
sqrt(abs(logFC.XY43.XX43))          1  35.32  35.321 1924.0158 < 2.2e-16 ***
  bias                                2   0.85   0.424   23.0809 9.663e-11 ***
  sqrt(abs(logFC.XY43.XX43)):bias     2   0.31   0.154    8.3743 0.0002314 ***
  Residuals                       25568 469.38   0.018                        
---
######

#scatter plot

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g43_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=g43_tau, xName='abslogFC.XY43.XX43',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC.XY46.XX46)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_g43_sb_nocolors.pdf", width=8, height=8)
ggplot2.scatterplot(data=g43_tau, xName='abslogFC.XY43.XX43',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="LogFC.XY43.XX43)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

########
g43_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/G43/Amg43_dnds_tau_match_sorted_fi.txt", header = TRUE)
str(g43_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_amm_g43_dnds.pdf", width=8, height=8)
ggplot(g43_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

#wilcox test for dNdS
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias=='male'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='female'],exact = FALSE) 
#W = 66, p-value = 0.9094
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias=='male'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 14262, p-value = 0.2533  ## only 3 for male-biased genes
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias=='unbias'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='female'],exact = FALSE) 
#W = 86214, p-value = 6.255e-06

##############
#load data from gonads
##############

gonad_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/gonad/Amgonad_log1_sb_un_tau_sorted_fi.txt", header = TRUE)
str(gonad_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_gonad_sb.pdf", width=8, height=8)
ggplot(gonad_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0, 1) +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

y2 <- lm(sqrt(tau)~sqrt(abs(logFC.XYtestis.XXovary))*bias, data=gonad_tau)
anova(y2)

####
Response: sqrt(tau)
Df Sum Sq Mean Sq F value    Pr(>F)    
sqrt(abs(logFC.XYtestis.XXovary))          1  84.83  84.833 5116.84 < 2.2e-16 ***
  bias                                       2  20.99  10.493  632.91 < 2.2e-16 ***
  sqrt(abs(logFC.XYtestis.XXovary)):bias     2   0.91   0.456   27.50 1.182e-12 ***
  Residuals                              21111 350.00   0.017                      
---
###

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_gonad_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau, xName='abslogFC.XYtestis.XXovary',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC(XYtestis/XXovary)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_gonad_sb.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau, xName='abslogFC.XYtestis.XXovary',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="LogFC(XYtestis/XXovary)", y="Tau") +
  scale_fill_manual(values = c("grey40"))

dev.off()

wilcox.test(gonad_tau$tau[gonad_tau$bias=='male'],gonad_tau$tau[gonad_tau$bias=='unbias'],exact = FALSE) 
#W = 33218000, p-value < 2.2e-16
wilcox.test(gonad_tau$tau[gonad_tau$bias=='female'],gonad_tau$tau[gonad_tau$bias=='unbias'],exact = FALSE) 
#W = 38308000, p-value < 2.2e-16
wilcox.test(gonad_tau$tau[gonad_tau$bias=='male'],gonad_tau$tau[gonad_tau$bias=='female'],exact = FALSE) 
#W = 10568000, p-value < 2.2e-16

##########
#####gonad tau and dnds
##########
gonad_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/gonad/Amgonad_log1_sb_un_tau_dnds_match_fi.txt", header = TRUE)
str(gonad_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newgonad_sb_dnds.pdf", width=8, height=8)
ggplot(gonad_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,1) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###linear model to investigate whether evolutionary rate is explained by tau or sex bias, or the interaction.
y4 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, data=gonad_tau_dnds)
anova(y4)

#######
Df Sum Sq Mean Sq  F value    Pr(>F)    
sqrt(tau)         1  1.191 1.19103 128.7234 < 2.2e-16 ***
  bias              2  0.002 0.00088   0.0949  0.909449    
sqrt(tau):bias    2  0.089 0.04437   4.7956  0.008299 ** 
  Residuals      5774 53.425 0.00925                       
##########

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 2430100, p-value = 0.002517

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1855500, p-value = 0.06863

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],exact = FALSE) 
#W = 1349300, p-value = 0.3569

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_dnds_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_dnds.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="dNdS", y="Tau") +
  scale_fill_manual(values = c("grey40"))

dev.off()

cor.test(gonad_tau_dnds$dNdS, gonad_tau_dnds$tau, method=c("pearson"))

######
data:  gonad_tau_dnds$dNdS and gonad_tau_dnds$tau
t = 9.4761, df = 5778, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.09823822 0.14901178
sample estimates:
  cor 
0.123706 
##########

##############
#load data from folder Brain tissues.
##############
brain_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/brain/Ambrain_log1_sb_un_tau_fi.txt", header = TRUE)
str(brain_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_amm_brain.pdf", width=8, height=8)
ggplot(brain_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

###some stats
wilcox.test(brain_tau$tau[brain_tau$bias=='female'],brain_tau$tau[brain_tau$bias=='unbias'],exact = FALSE) 
#W = 1642100, p-value < 2.2e-16

wilcox.test(brain_tau$tau[brain_tau$bias=='male'],brain_tau$tau[brain_tau$bias=='unbias'],exact = FALSE) 
#W = 1328800, p-value = 2.902e-07

wilcox.test(brain_tau$tau[brain_tau$bias=='male'],brain_tau$tau[brain_tau$bias=='female'],exact = FALSE) 
#W = 1296, p-value = 3.262e-08

###dNdS and tau in Brain tissues
brain_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/brain/Ambrain_tau_dnds_match_fi.txt", header = TRUE)
str(brain_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/brain_newsbun_dnds.pdf", width=8, height=8)
ggplot(brain_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,1) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###some stats
wilcox.test(brain_tau_dnds$dNdS[brain_tau_dnds$bias=='female'],brain_tau_dnds$dNdS[brain_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 24332, p-value = 0.1197
wilcox.test(brain_tau_dnds$dNdS[brain_tau_dnds$bias=='male'],brain_tau_dnds$dNdS[brain_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 103310, p-value = 0.1372
wilcox.test(brain_tau_dnds$dNdS[brain_tau_dnds$bias=='female'],brain_tau_dnds$dNdS[brain_tau_dnds$bias=='male'],exact = FALSE) 
#W = 74, p-value = 0.04997

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_dnds_brain_colors.pdf", width=8, height=8)
ggplot2.scatterplot(data=brain_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS", y="Tau", color ="bias")
dev.off()


########
qqnorm(sqrt(abs(brain_tau$tau)))
hist(sqrt(abs(brain_tau$tau)))
qqnorm(sqrt(abs(brain_tau_dnds$dNdS)))
hist(sqrt(abs(brain_tau_dnds$dNdS)))

y2 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, brain_tau_dnds)
anova(y2)

##########
Df Sum Sq Mean Sq F value Pr(>F)    
sqrt(tau)         1  0.645 0.64482 70.8861 <2e-16 ***
  bias              2  0.034 0.01724  1.8947 0.1504    
sqrt(tau):bias    2  0.034 0.01711  1.8808 0.1525    
Residuals      6830 62.129 0.00910  
##########

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_brain_exp.pdf", width=8, height=8)
ggplot2.scatterplot(data=brain_tau_dnds, xName='abslogFC.XYbrain.Xxrain',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="abs(logFC.XYbrain.XXrain)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

##############
#load data from folder Liver tissues.
##############
liver_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/liver/Amliver_log1_sb_un_tau_sorted_fi.txt", header = TRUE)
str(liver_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/newtau_amm_liver.pdf", width=8, height=8)
ggplot(liver_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

###some stats
wilcox.test(liver_tau$tau[liver_tau$bias=='female'],liver_tau$tau[liver_tau$bias=='unbias'],exact = FALSE) 
#W = 1890600, p-value < 2.2e-16
wilcox.test(liver_tau$tau[liver_tau$bias=='male'],liver_tau$tau[liver_tau$bias=='unbias'],exact = FALSE) 
#W = 1286200, p-value = 7.715e-11
wilcox.test(liver_tau$tau[liver_tau$bias=='female'],liver_tau$tau[liver_tau$bias=='male'],exact = FALSE) 
#W = 7687, p-value = 0.1625

####
liver_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/liver/Amliver_tau_dnds_match_fi.txt", header = TRUE)
str(liver_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/liver_newsbun_dnds.pdf", width=8, height=8)
ggplot(liver_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,1) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###some stats
wilcox.test(liver_tau_dnds$dNdS[liver_tau_dnds$bias=='female'],liver_tau_dnds$dNdS[liver_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 81166, p-value = 0.0796
wilcox.test(liver_tau_dnds$dNdS[liver_tau_dnds$bias=='male'],liver_tau_dnds$dNdS[liver_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 123380, p-value = 0.005197
wilcox.test(liver_tau_dnds$dNdS[liver_tau_dnds$bias=='female'],liver_tau_dnds$dNdS[liver_tau_dnds$bias=='male'],exact = FALSE) 
#W = 426.5, p-value = 0.7359

####
qqnorm(sqrt(abs(liver_tau$tau)))
hist(sqrt(abs(liver_tau$tau)))
qqnorm(sqrt(abs(liver_tau_dnds$dNdS)))
hist(sqrt(abs(liver_tau_dnds$dNdS)))

y3 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, liver_tau_dnds)
anova(y3)

##
Df Sum Sq Mean Sq  F value  Pr(>F)    
sqrt(tau)         1  1.910 1.91017 212.7090 < 2e-16 ***
  bias              2  0.082 0.04080   4.5432 0.01068 *  
  sqrt(tau):bias    2  0.051 0.02525   2.8121 0.06016 .  
Residuals      5453 48.969 0.00898                    
---
  ##
  
  #with color
  pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_dnds_liver_colors.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS", y="Tau", color ="bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_exp_liver_colors.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='abslogFC.XYliver.Xxliver',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="abs(logFC.XYliver.XXliver)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_liver_dnds.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="dNdS", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/new_tau/scatter_abs_newtau_exp_liver_nocolors.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='abslogFC.XYliver.Xxliver',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="abs(logFC.XYliver.XXliver)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()
