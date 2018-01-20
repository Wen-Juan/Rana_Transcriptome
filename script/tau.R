#intall R packages and load libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
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

##############
#load data from folder G46
##############

g46_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_tau/G46_sb_un_tau_match_fi.txt", header = TRUE)
str(g46_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g46.pdf", width=8, height=8)

ggplot(g46_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
geom_boxplot() +
labs(x='Sex bias', y='Tau')

dev.off()

####liner model to investigate whether tau is correlated with gene expression and sex bias.
g46_tau_sub <- subset(g46_tau,g46_tau$bias!='unbias')
str(g46_tau_sub)
qqnorm(sqrt(abs(g46_tau_sub$logFC.XY46.XX46)))
hist(sqrt(abs(g46_tau_sub$logFC.XY46.XX46)))
y <- lm(tau~sqrt(abs(g46_tau_sub$logFC.XY46.XX46))*bias, g46_tau_sub)
anova(y)

##########
##############without unbias genes
##########
#Response: tau
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#sqrt(abs(g46_tau_sub$logFC.XY46.XX46))         1 54.872  54.872 3407.645 < 2.2e-16 ***
#  bias                                           1  1.530   1.530   94.997 < 2.2e-16 ***
#  sqrt(abs(g46_tau_sub$logFC.XY46.XX46)):bias    1  2.275   2.275  141.265 < 2.2e-16 ***
#  Residuals                                   4674 75.264   0.016 

#############
#############with unbias genes
#Response: tau
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#sqrt(abs(g46_tau$logFC.XY46.XX46))          1 108.38 108.379 3855.423 < 2.2e-16 ***
#  bias                                        2   5.12   2.562   91.155 < 2.2e-16 ***
#  sqrt(abs(g46_tau$logFC.XY46.XX46)):bias     2   6.91   3.454  122.879 < 2.2e-16 ***
#  Residuals                               24400 685.90   0.028  
############
###########

##investigate whether dnds is influences by tau or sex bias, or the interaction between the two.
g46_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_tau/G46_sb_un_tau_dnds.txt", header = TRUE)
str(g46_tau_dnds)

hist(g46_tau_dnds$tau) #look quite normal distribution
y1 <- lm(dNdS~tau*bias, g46_tau_dnds)
anova(y1)

###########
#Response: dNdS
#Df  Sum Sq Mean Sq  F value    Pr(>F)    
#  tau          1  0.4084 0.40841 119.9465 < 2.2e-16 ***
#  bias         2  0.0652 0.03258   9.5674 7.095e-05 ***
#  tau:bias     2  0.1277 0.06386  18.7562 7.540e-09 ***
#  Residuals 6589 22.4353 0.00340  
###########

y2 <- glm(dNdS~tau*bias, family=quasi, g46_tau_dnds)
summary(y2)
#######
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)     0.038071   0.004829   7.884 3.69e-15 ***
  tau             0.092544   0.009279   9.974  < 2e-16 ***
  biasmale        0.037546   0.018415   2.039   0.0415 *  
  biasunbias      0.031121   0.005254   5.923 3.33e-09 ***
  tau:biasmale   -0.033100   0.033230  -0.996   0.3192    
tau:biasunbias -0.062678   0.010277  -6.099 1.13e-09 ***
#########  

#scatter plot
abs_g46_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_tau/G46_abssb_un_tau_match_fi.txt", header = TRUE)
str(abs_g46_tau)

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g46_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC.XY46.XX46)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g46_sb.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="LogFC.XY46.XX46)", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

#stats
wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='unbias'],exact = FALSE) 
#W = 9496800, p-value < 2.2e-16

wilcox.test(g46_tau$tau[g46_tau$bias=='male'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 1484700, p-value = 0.03636

wilcox.test(g46_tau$tau[g46_tau$bias=='unbias'],g46_tau$tau[g46_tau$bias=='female'],exact = FALSE)
#W = 26901000, p-value < 2.2e-16

###########
#for sex-biased genes and tau in gonad tissues
##########
##boxplot
gonad_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Gonad_tau/gonad_fc1_sb_unbias_tau.txt", header = TRUE)
str(gonad_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_gonad_sb.pdf", width=8, height=8)
ggplot(gonad_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  labs(x='Sex bias', y='Tau') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

dev.off()

y <- lm(tau~sqrt(abs(logFC.XYtestis.XXovary))*bias, gonad_tau)
anova(y)
##
##
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#sqrt(abs(logFC.XYtestis.XXovary))          1 165.94 165.945 6088.711 < 2.2e-16 ***
#  bias                                       2  34.67  17.333  635.955 < 2.2e-16 ***
#  sqrt(abs(logFC.XYtestis.XXovary)):bias     2   4.06   2.028   74.402 < 2.2e-16 ***
#  Residuals                              21111 575.37   0.027  
##


##correlation scatter plot
abs_gonad_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Gonad_tau/gonad_fc1_abssb_unbias_tau.txt", header = TRUE)
str(abs_gonad_tau)

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_gonad_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_gonad_tau, xName='abslogFC.XYtestis.XXovary',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC(XYtestis/XXovary)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_gonad_sb.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_gonad_tau, xName='abslogFC.XYtestis.XXovary',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="LogFC(XYtestis/XXovary)", y="Tau") +
  scale_fill_manual(values = c("grey40"))

 dev.off()

wilcox.test(gonad_tau$tau[gonad_tau$bias=='male'],gonad_tau$tau[gonad_tau$bias=='unbias'],exact = FALSE) 
#W = 33823000, p-value < 2.2e-16
wilcox.test(gonad_tau$tau[gonad_tau$bias=='female'],gonad_tau$tau[gonad_tau$bias=='unbias'],exact = FALSE) 
#W = 38865000, p-value < 2.2e-16
wilcox.test(gonad_tau$tau[gonad_tau$bias=='male'],gonad_tau$tau[gonad_tau$bias=='female'],exact = FALSE) 
#W = 10587000, p-value < 2.2e-16

##########
#####gonad tau and dnds
##########
gonad_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Gonad_tau/gonad_fc1_tau_dnds_subset_sorted.txt", header = TRUE)
str(gonad_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/gonad_sb_dnds.pdf", width=8, height=8)
ggplot(gonad_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,0.5) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###linear model to investigate whether evolutionary rate is explained by tau or sex bias, or the interaction.
y1 <- lm(dNdS ~ sqrt(tau) * bias, gonad_tau_dnds)
anova(y1)
#######
#Response: dNdS
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#sqrt(tau)         1  0.2805 0.280491 80.0801 < 2.2e-16 ***
#  bias              2  0.0020 0.001013  0.2892  0.748840    
#sqrt(tau):bias    2  0.0378 0.018880  5.3904  0.004583 ** 
#  Residuals      5774 20.2242 0.003503 
########

y1a <- glm(dNdS ~ tau * bias,family=quasipoisson, gonad_tau_dnds)
summary(y1a)

#####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)    -2.77772    0.05513 -50.389  < 2e-16 ***
  tau             0.62930    0.10043   6.266 3.97e-10 ***
  biasmale        0.22170    0.07373   3.007  0.00265 ** 
  biasunbias      0.09160    0.06541   1.401  0.16141    
tau:biasmale   -0.48594    0.13841  -3.511  0.00045 ***
  tau:biasunbias -0.21891    0.12770  -1.714  0.08652 . 
#####

y2 <- glm(dNdS ~ tau * logFC.XYtestis.XXovary, family=quasipoisson, gonad_tau_dnds)
summary(y2)

#########
Coefficients:
  Estimate Std. Error  t value Pr(>|t|)    
(Intercept)                -2.67217    0.02510 -106.456  < 2e-16 ***
  tau                         0.38679    0.05060    7.644 2.45e-14 ***
  logFC.XYtestis.XXovary      0.04836    0.01376    3.514 0.000445 ***
  tau:logFC.XYtestis.XXovary -0.09222    0.02108   -4.374 1.24e-05 ***
#########  


wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 2430100, p-value = 0.002517

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1855500, p-value = 0.06863

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],exact = FALSE) 
#W = 1349300, p-value = 0.3569

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_dnds_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_dnds.pdf", width=8, height=8)
ggplot2.scatterplot(data=gonad_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="dNdS", y="Tau") +
  scale_fill_manual(values = c("grey40"))

dev.off()

cor.test(gonad_tau_dnds$dNdS, gonad_tau_dnds$tau, method=c("pearson"))

######
#t = 8.0396, df = 5778, p-value = 1.085e-15
# cor 
#0.1051791
######

cor.test(gonad_tau_dnds$abslogFC.XYtestis.XXovary, gonad_tau_dnds$tau, method=c("pearson"))
#t = 35.609, df = 5778, p-value < 2.2e-16
#0.4242161

cor.test(gonad_tau_dnds$abslogFC.XYtestis.XXovary, gonad_tau_dnds$dNdS, method=c("pearson"))
#t = 8.0656, df = 5778, p-value = 8.792e-16
#     cor 
#0.1055151 

cor.test(gonad_tau_dnds$abslogFC.XYtestis.XXovary, gonad_tau_dnds$dNdS, method=c("pearson"))