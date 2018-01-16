#intall R packages and load libraries
install.packages("ggplot2")

library("ggplot2")


#set working directory
datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/'
results <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/'

#loading data
Amm_all_TPM <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Amm_all_TPM.txt",header = T)
str(Amm_all_TPM)

###tissue number
######developmental stages, 7 tissues in total: G23, G27, G31, G43_female(G43_XX), G43_male(G43_XY), G46_female(G46_XX), G46_male(G46_XY) 
######adult tissues, 4 tissues in total: testis, ovary, brain, liver

###function requires data frame with expression values
#mean values between replicates are calculated.

Amm_all_TPM$G23_mean <- (Amm_all_TPM$G23 + Amm_all_TPM$G23.1 + Amm_all_TPM$G23.2 +Amm_all_TPM$G23.3 +Amm_all_TPM$G23.4 +Amm_all_TPM$G23.5)/6
Amm_all_TPM$G27_mean <- (Amm_all_TPM$G27 + Amm_all_TPM$G27.1 + Amm_all_TPM$G27.2 +Amm_all_TPM$G27.3 +Amm_all_TPM$G27.4 +Amm_all_TPM$G27.5)/6
Amm_all_TPM$G31_mean <- (Amm_all_TPM$G31 + Amm_all_TPM$G31.1 + Amm_all_TPM$G31.2 +Amm_all_TPM$G31.3 +Amm_all_TPM$G31.4 +Amm_all_TPM$G31.5)/6
Amm_all_TPM$G43_M_mean <- (Amm_all_TPM$G43_M + Amm_all_TPM$G43_M.1 + Amm_all_TPM$G43_M.2)/3
Amm_all_TPM$G43_F_mean <- (Amm_all_TPM$G43_F + Amm_all_TPM$G43_F.1 + Amm_all_TPM$G43_F.2)/3
Amm_all_TPM$G46_F_mean <- (Amm_all_TPM$G46_F + Amm_all_TPM$G46_F.1)/2
Amm_all_TPM$G46_M_mean <- (Amm_all_TPM$G46_M + Amm_all_TPM$G46_M.1 + Amm_all_TPM$G46_M.2)/3
Amm_all_TPM$Testis_mean <- (Amm_all_TPM$Testis + Amm_all_TPM$Testis.1 + Amm_all_TPM$Testis.2 +Amm_all_TPM$Testis.3 +Amm_all_TPM$Testis.4)/5
Amm_all_TPM$Ovary_mean <- (Amm_all_TPM$Ovary + Amm_all_TPM$Ovary.1 + Amm_all_TPM$Ovary.2 +Amm_all_TPM$Ovary.3 +Amm_all_TPM$Ovary.4)/5
Amm_all_TPM$Brain_mean <- (Amm_all_TPM$Brain + Amm_all_TPM$Brain.1 + Amm_all_TPM$Brain.2 +Amm_all_TPM$Brain.3 +Amm_all_TPM$Brain.4 +Amm_all_TPM$Brain.5+Amm_all_TPM$Brain.6+Amm_all_TPM$Brain.7+Amm_all_TPM$Brain.8+Amm_all_TPM$Brain.9)/10
Amm_all_TPM$Liver_mean <- (Amm_all_TPM$Liver + Amm_all_TPM$Liver.1 + Amm_all_TPM$Liver.2 +Amm_all_TPM$Liver.3 +Amm_all_TPM$Liver.4 +Amm_all_TPM$Liver.5+Amm_all_TPM$Liver.6+Amm_all_TPM$Liver.7+Amm_all_TPM$Liver.8)/9
Amm_all_TPM$all_mean <- (Amm_all_TPM$G23_mean+Amm_all_TPM$G27_mean+Amm_all_TPM$G31_mean+Amm_all_TPM$G43_M_mean+Amm_all_TPM$G43_F_mean+Amm_all_TPM$G46_F_mean+Amm_all_TPM$G46_M_mean+Amm_all_TPM$Testis_mean+Amm_all_TPM$Ovary_mean+Amm_all_TPM$Liver_mean+Amm_all_TPM$Brain_mean)/11


#select TPM>=1 in at least one tissue out of 11 tissues, average TPM>=0.1
Amm_all_TPM <- Amm_all_TPM[Amm_all_TPM$all_mean>=0.1,]
a1 <- Amm_all_TPM
m <- 10
rownames(a1) <- a1[,1]
a1 <-a1[,-1]
a1 <- round(a1,3)
dim1 <- dim(a1)[1]
a2 <- a1[,c(62,63:72)]
dim1 <- dim(a2)[1]

s <-matrix(0, nrow=dim(a2)[1], ncol=2, dimnames = list(rownames(a2),c("Tissue_specific_index", "Tissue_name")))

#with log transformation, reduce effects of extreme data and variation inbetween.
for(i in 1:dim1){ 
  s[i,1] <- (sum(1 - log(1+(a2[i,]))/max(log(1+(a2[i,])))))/ m   
  s[i,2] <- names(which.max(log(1+(a2[i,]))))
}

write.table(s, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm.txt",row.names=TRUE, sep="\t", quote= FALSE)

#plot tau
tau_amm <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm.txt", header=T)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_logtransf.pdf", width=8, height=8)
ggplot(tau_amm, aes(Tissue_specific_index)) +
  geom_histogram(binwidth = 0.01,fill=I("blue"),col=I("blue"))
dev.off()

#without log tansformation, data might be with higher variation
for(i in 1:dim1){ 
  s[i,1] <- (sum(1 - (a2[i,])/max(a2[i,])))/ m   
  s[i,2] <- names(which.max(a2[i,]))
}

write.table(s, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_withoutlogtransf.txt",row.names=TRUE, sep="\t", quote= FALSE)

#make ggplot for non-log transformed tau distribution
tau_amm_nolog <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_withoutlogtransf.txt", header=T)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_withoutlogtransf.pdf", width=8, height=8)
ggplot(tau_amm_nolog, aes(Tissue_specific_index)) +
  geom_histogram(binwidth = 0.01,fill=I("blue"),col=I("blue"))
dev.off()

#subset the tissue specific genes (high tau values, >=0.9), and also genes with more pleotropic effects (tau<0.2).
s2 <- subset(tau_amm,tau_amm[,2] >= 0.9)
write.table(s2, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_log_TSI0.9.txt",row.names=TRUE, sep="\t", quote= FALSE)
s3 <- subset(tau_amm,tau_amm[,2]<=0.2)
write.table(s3, file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_log_TSI0.2.txt",row.names=TRUE, sep="\t", quote= FALSE)

#load data from folder
g46_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_sb_un_tau_match_fi.txt", header = TRUE)
str(g46_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g46.pdf", width=8, height=8)

ggplot(g46_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
geom_boxplot() +
labs(x='Sex bias', y='Tau')

dev.off()

#stats
wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='unbias'],exact = FALSE) 
#W = 9496800, p-value < 2.2e-16

wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 1484700, p-value = 0.03636

wilcox.test(g46_tau$tau[g46_tau$bias=='unbias'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 26901000, p-value < 2.2e-16
