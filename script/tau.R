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
a2 <- a1[,c(1,63:73)]
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
ylim(0,1.25) +
scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

########
cor.test(sqrt(abs(g46_tau$logFC.XY46.XX46)),sqrt(g46_tau$tau),
         method = "spearman",
         conf.level = 0.95)
####
data:  sqrt(abs(g46_tau$logFC.XY46.XX46)) and sqrt(g46_tau$tau)
S = 1.7071e+12, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2954389 
####

####liner model to investigate whether tau is correlated with gene expression and sex bias.
qqnorm(sqrt(g46_tau$tau))
hist(sqrt(g46_tau$tau))

y <- lm(sqrt(tau)~g46_tau$logFC.XY46.XX46*bias, data=g46_tau)
summary(y)
anova(y)  #does not provide p values for some reason

##########
Df Sum Sq Mean Sq F value    Pr(>F)    
g46_tau$logFC.XY46.XX46          1  18.36 18.3641 1016.32 < 2.2e-16 ***
  bias                             2  11.62  5.8113  321.62 < 2.2e-16 ***
  g46_tau$logFC.XY46.XX46:bias     2  26.25 13.1267  726.47 < 2.2e-16 ***
  Residuals                    24400 440.89  0.0181   

y1 <- glm(tau~sqrt(abs(g46_tau$logFC.XY46.XX46))*bias, family=binomial, data=g46_tau)
summary(y1)
anova(y) 
#######

##investigate whether dnds is influences by tau or sex bias, or the interaction between the two.
###
g46_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_tau/G46_sb_shareunbias_tau_dnds.txt", header = TRUE)
str(g46_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g46_dnds_shareunbias.pdf", width=8, height=8)
ggplot(g46_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='male'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 145600, p-value = 1.843e-09
wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='female'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1362300, p-value = 0.00253
wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='male'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='female'],exact = FALSE) 
#W = 152450, p-value = 3.095e-06
g46sub_tau_dnds <- subset(g46_tau_dnds, g46_tau_dnds$bias!='unbias')
str(g46sub_tau_dnds)

###tau
ggplot(g46_tau_dnds, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1.25) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))

wilcox.test(g46_tau_dnds$tau[g46_tau_dnds$bias=='male'],g46_tau_dnds$tau[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 184270, p-value < 2.2e-16
wilcox.test(g46_tau_dnds$tau[g46_tau_dnds$bias=='female'],g46_tau_dnds$tau[g46_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1991900, p-value < 2.2e-16
wilcox.test(g46_tau_dnds$dNdS[g46_tau_dnds$bias=='male'],g46_tau_dnds$dNdS[g46_tau_dnds$bias=='female'],exact = FALSE) 
#W = 152450, p-value = 3.095e-06

qqnorm(sqrt(g46_tau_dnds$tau)[1:5000]) #not normal distribution
shapiro.test(log(g46_tau_dnds$tau[1:5000])) #not normal distribution

y1 <- lm(dNdS~sqrt(tau)*bias, g46_tau_dnds)
anova(y1)

###########
#Response: dNdS
#Df  Sum Sq Mean Sq  F value    Pr(>F)    
#  sqrt(tau)         1  0.4781 0.47813 140.8440 < 2.2e-16 ***
#  bias              2  0.0643 0.03216   9.4725 7.799e-05 ***
#  sqrt(tau):bias    2  0.1263 0.06317  18.6083 8.735e-09 ***
#  Residuals      6589 22.3678 0.00339  
###########

################
#########if removing unbiased genes,
###############
y2 <- lm(dNdS~sqrt(tau)*bias, g46sub_tau_dnds)
anova(y2)

############
#Response: dNdS
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#  sqrt(tau)         1 0.3712 0.37124 115.7551 < 2.2e-16 ***
#  bias              1 0.0534 0.05341  16.6526 4.682e-05 ***
#  sqrt(tau):bias    1 0.0034 0.00336   1.0475    0.3062    
#   Residuals      1827 5.8594 0.00321 
############

qqnorm(sqrt(g46_tau_dnds$dNdS))

y1a <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, g46_tau_dnds)
anova(y1a)
################
Analysis of Variance Table

Response: sqrt(dNdS)
Df Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)         1  1.607 1.60698 175.214 < 2.2e-16 ***
  bias              2  0.197 0.09855  10.745 2.192e-05 ***
  sqrt(tau):bias    2  0.403 0.20170  21.992 3.024e-10 ***
  Residuals      6589 60.431 0.00917                      
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#################
summary(y1a)
###############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)           0.10758    0.01506   7.143 1.01e-12 ***
  sqrt(tau)             0.23680    0.02136  11.088  < 2e-16 ***
  biasmale              0.10490    0.05529   1.897   0.0578 .  
biasunbias            0.10901    0.01620   6.730 1.84e-11 ***
  sqrt(tau):biasmale   -0.09704    0.07573  -1.282   0.2001    
sqrt(tau):biasunbias -0.15363    0.02322  -6.617 3.96e-11 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############

##if removing unbiased genes
g46_tau_dnds_sub <- subset (g46_tau_dnds,g46_tau_dnds$bias!='unbias')
str(g46_tau_dnds_sub)

y2a <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, g46_tau_dnds_sub)
anova(y2a)
###############
Response: sqrt(dNdS)
Df  Sum Sq Mean Sq  F value    Pr(>F)    
sqrt(tau)         1  1.2116 1.21159 132.5928 < 2.2e-16 ***
  bias              1  0.1634 0.16337  17.8782 2.471e-05 ***
  sqrt(tau):bias    1  0.0151 0.01506   1.6484    0.1993    
Residuals      1827 16.6946 0.00914   
#############
summary(y2a)
#############
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.10758    0.01503   7.156  1.2e-12 ***
  sqrt(tau)           0.23680    0.02132  11.109  < 2e-16 ***
  biasmale            0.10490    0.05518   1.901   0.0575 .  
sqrt(tau):biasmale -0.09704    0.07559  -1.284   0.1993    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#############

#scatter plot
abs_g46_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G46_tau/G46_abssb_un_tau_match_fi.txt", header = TRUE)
str(abs_g46_tau)

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g46_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC.XY46.XX46", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g46_sb.pdf", width=8, height=8)
p2 <- ggplot2.scatterplot(data=abs_g46_tau, xName='abslogFC.XY46.XX46',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="G46: Log2(male/female)", y="Tissue specificity") +
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
  ylim(0, 1) +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

gonad_tau_sub <- subset(gonad_tau,gonad_tau$chr=='Chr01')
shapiro.test(gonad_tau_sub$tau) #not significant

y2 <- glm(sqrt(tau)~sqrt(abs(logFC.XYtestis.XXovary))*bias, family = binomial, data=gonad_tau)
summary(y2)

##
##
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                  -0.45131    0.11497  -3.926 8.65e-05 ***
#  sqrt(abs(logFC.XYtestis.XXovary))             0.91428    0.07619  11.999  < 2e-16 ***
#  biasmale                                     -0.09787    0.16565  -0.591    0.555    
# biasunbias                                    0.53591    0.12671   4.229 2.34e-05 ***
#  sqrt(abs(logFC.XYtestis.XXovary)):biasmale   -0.15882    0.10971  -1.448    0.148    
# sqrt(abs(logFC.XYtestis.XXovary)):biasunbias -0.43344    0.10830  -4.002 6.28e-05 ***
##

##############
#load data from folder G43
##############
g43_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G43_tau/Am43_log1_sb_sorted.txt", header = TRUE)
str(g43_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g43.pdf", width=8, height=8)
ggplot(g43_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

#wilcox test for tau
wilcox.test(g43_tau$tau[g43_tau$bias=='male'],g43_tau$tau[g43_tau$bias=='female'],exact = FALSE) 
#W = 664, p-value = 5.495e-10
wilcox.test(g43_tau$tau[g43_tau$bias=='unbias'],g43_tau$tau[g43_tau$bias=='female'],exact = FALSE) 
#W = 191360, p-value < 2.2e-16
wilcox.test(g43_tau$tau[g43_tau$bias=='unbias'],g43_tau$tau[g43_tau$bias=='male'],exact = FALSE) 
#W = 186350, p-value = 1.289e-06

####some stats

str(g43_tau)
qqnorm(sqrt(abs(g43_tau$logFC.XY43.XX43)))
hist(sqrt(abs(g43_tau$logFC.XY43.XX43)))
qqnorm(sqrt(abs(g43_tau$tau)))
hist(sqrt(abs(g43_tau$tau)))

cor.test(sqrt(g43_tau$abslogFC.XY43.XX43),sqrt(g43_tau$tau),
         method = "spearman",
         conf.level = 0.95)

########
data:  sqrt(g43_tau$abslogFC.XY43.XX43) and sqrt(g43_tau$tau)
S = 2.1933e+12, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.2131365 
#######

y <- lm(sqrt(tau)~sqrt(abs(logFC.XY43.XX43))*bias, data=g43_tau)
anova(y)

######
#> anova(y)
#Analysis of Variance Table
#Response: sqrt(tau)
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#sqrt(abs(logFC.XY43.XX43))          1  35.31  35.307 1896.228 < 2.2e-16 ***
#  bias                                2   0.60   0.300   16.094 1.035e-07 ***
#  sqrt(abs(logFC.XY43.XX43)):bias     2   0.38   0.189   10.137 3.974e-05 ***
#  Residuals                       25567 476.04   0.019                       
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
######

#scatter plot ###
####

absg43tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G43_tau/Am43_log1_sb_sorted.txt", header = TRUE)
str(absg43tau)

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g43_sb_colors1.pdf", width=8, height=8)
ggplot2.scatterplot(data=abs_g43_tau, xName='abslogFC.XY43.XX43',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="LogFC.XY46.XX46)", y="Tau", color ="bias")
dev.off()


#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_g43_sb_nocolors.pdf", width=8, height=8)
p1 <- ggplot2.scatterplot(data=absg43tau, xName='abslogFC.XY43.XX43',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="G43: Log2(male/female)", y="Tissue specificity") +
  scale_fill_manual(values = c("grey40"))
dev.off()

########
g43_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/G43_tau/Am43_log1_sb_shareunbias_dnds_sorted.txt", header = TRUE)
str(g43_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g43_dnds_shareunbias.pdf", width=8, height=8)
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
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias!='unbias'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 49395, p-value = 4.594e-07
#W = 14262, p-value = 0.2533  ## only 3 for male-biased genes
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias=='unbias'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='female'],exact = FALSE) 
#W = 17852, p-value = 9.837e-07
wilcox.test(g43_tau_dnds$dNdS[g43_tau_dnds$bias=='unbias'],g43_tau_dnds$dNdS[g43_tau_dnds$bias=='male'],exact = FALSE) 
#W = 1333, p-value = 0.2119

###tau
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_g43_tau_shareunbias.pdf", width=8, height=8)

ggplot(g43_tau_dnds, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

wilcox.test(g43_tau_dnds$tau[g43_tau_dnds$bias=='male'],g43_tau_dnds$tau[g43_tau_dnds$bias=='female'],exact = FALSE) 
#W = 4, p-value = 0.007771
wilcox.test(g43_tau_dnds$tau[g43_tau_dnds$bias!='unbias'],g43_tau_dnds$tau[g43_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 67275, p-value < 2.2e-16
#W = 14262, p-value = 0.2533  ## only 3 for male-biased genes
wilcox.test(g43_tau_dnds$tau[g43_tau_dnds$bias=='unbias'],g43_tau_dnds$tau[g43_tau_dnds$bias=='female'],exact = FALSE) 
#W = 633, p-value < 2.2e-16
wilcox.test(g43_tau_dnds$tau[g43_tau_dnds$bias=='unbias'],g43_tau_dnds$tau[g43_tau_dnds$bias=='male'],exact = FALSE) 
#W = 672, p-value = 0.03446


y <- lm(sqrt(dNdS)~sqrt(tau)*bias, data=g43_tau_dnds)
anova(y)

###
Response: sqrt(dNdS)
Df  Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)         1  0.5102 0.51021  55.685 1.405e-13 ***
  bias              2  0.1027 0.05136   5.606   0.00375 ** 
  sqrt(tau):bias    2  0.0535 0.02677   2.922   0.05412 .  
Residuals      1563 14.3208 0.00916    
###


y1 <- lm(sqrt(dNdS)~sqrt(tau) + bias, data=g43_tau_dnds)
anova(y1)

###
Response: sqrt(dNdS)
Df  Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)    1  0.5102 0.51021 55.5491 1.501e-13 ***
  bias         2  0.1027 0.05136  5.5923  0.003801 ** 
  Residuals 1565 14.3743 0.00918   
###

##############
#load data from gonads
##############

###If removing unbiased genes, 
gonad_tau_sub <- subset(gonad_tau,gonad_tau$bias!='unbias')

hist(sqrt(abs_gonad_tau$abslogFC.XYtestis.XXovary))
hist(sqrt(abs_gonad_tau$tau))

cor.test(sqrt(abs_gonad_tau$abslogFC.XYtestis.XXovary),sqrt(abs_gonad_tau$tau),
         method = "spearman",
         conf.level = 0.95)

#######
data:  sqrt(abs_gonad_tau$abslogFC.XYtestis.XXovary) and sqrt(abs_gonad_tau$tau)
S = 9.5359e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.3924036 
#######

y3 <- lm(sqrt(tau)~sqrt(abs(logFC.XYtestis.XXovary))*bias, data=gonad_tau)
summary(y3)

#################
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                -0.45131    0.11497  -3.926 8.65e-05 ***
#  sqrt(abs(logFC.XYtestis.XXovary))           0.91428    0.07619  11.999  < 2e-16 ***
#  biasmale                                   -0.09787    0.16565  -0.591    0.555    
# sqrt(abs(logFC.XYtestis.XXovary)):biasmale -0.15882    0.10971  -1.448    0.148  
#################


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
p3 <- ggplot2.scatterplot(data=abs_gonad_tau, xName='abslogFC.XYtestis.XXovary',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="Gonad: Log2(male/female)", y="Tissue specificity") +
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
gonad_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Gonad_tau/gonad_fc1_tau_dnds_sb_shareunbias_sorted_fi.txt", header = TRUE)
str(gonad_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/gonad_sb_dnds_shareunbias.pdf", width=8, height=8)
ggplot(gonad_tau_dnds, aes(x=bias, y=dNdS, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,0.6) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###linear model to investigate whether evolutionary rate is explained by tau or sex bias, or the interaction.
y4 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, data=gonad_tau_dnds)
anova(y4)
##
Df Sum Sq Mean Sq  F value    Pr(>F)    
sqrt(tau)         1  1.191 1.19103 128.7234 < 2.2e-16 ***
  bias              2  0.002 0.00088   0.0949  0.909449    
sqrt(tau):bias    2  0.089 0.04437   4.7956  0.008299 ** 
  Residuals      5774 53.425 0.00925                       
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###

summary(y4)

#######
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)           0.18586    0.01392  13.348  < 2e-16 ***
  sqrt(tau)             0.12870    0.01954   6.585 4.94e-11 ***
  biasmale              0.05603    0.01805   3.104  0.00192 ** 
  biasunbias            0.02260    0.01618   1.397  0.16253    
sqrt(tau):biasmale   -0.08342    0.02576  -3.239  0.00121 ** 
  sqrt(tau):biasunbias -0.03435    0.02356  -1.458  0.14494    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##########

#removing unbiased genes.
gonad_tau_dnds_sub <- subset(gonad_tau_dnds,gonad_tau_dnds$bias!='unbias')

#check normal distribution  
qqnorm(sqrt(gonad_tau_dnds_sub$dNdS))
hist(sqrt(gonad_tau_dnds$dNdS))

#ad.test for large dataset
install.packages("nortest")
library(nortest)
ad.test(sqrt(gonad_tau_dnds_sub$dNdS))

#ks test for large dataset
ks.test(sqrt(gonad_tau_dnds_sub$dNdS), "pnorm")

y5 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, gonad_tau_dnds_sub)
summary(y5)
#####
Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.18586    0.01388  13.392  < 2e-16 ***
  sqrt(tau)           0.12870    0.01948   6.607 4.55e-11 ***
  biasmale            0.05603    0.01799   3.114  0.00186 ** 
  sqrt(tau):biasmale -0.08342    0.02567  -3.250  0.00117 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1   
#####

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1513800, p-value = 0.0001378

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1156200, p-value = 0.006921

wilcox.test(gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$dNdS[gonad_tau_dnds$bias=='male'],exact = FALSE) 
#W = 1349300, p-value = 0.3569

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/gonad_sb_tau_shareunbias.pdf", width=8, height=8)
ggplot(gonad_tau_dnds, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,1) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()


wilcox.test(gonad_tau_dnds$tau[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$tau[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 2207300, p-value < 2.2e-16

wilcox.test(gonad_tau_dnds$tau[gonad_tau_dnds$bias=='male'],gonad_tau_dnds$tau[gonad_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 1533000, p-value < 2.2e-16

wilcox.test(gonad_tau_dnds$tau[gonad_tau_dnds$bias=='female'],gonad_tau_dnds$tau[gonad_tau_dnds$bias=='male'],exact = FALSE) 
#W = 1495500, p-value = 2.101e-10

###########
#remove unbiased genes.
##########
gonadsub_tau_dnds <- subset(gonad_tau_dnds, gonad_tau_dnds$bias!='unbias')
str(gonadsub_tau_dnds)
y2 <- lm(dNdS ~ sqrt(tau) * bias, gonadsub_tau_dnds)
anova(y2)

#Response: dNdS
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#  sqrt(tau)         1  0.1280 0.127956 35.2817 3.152e-09 ***
#  bias              1  0.0014 0.001430  0.3944  0.530034    
#  sqrt(tau):bias    1  0.0376 0.037630 10.3757  0.001289 ** 
#  Residuals      3277 11.8847 0.003627    


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


##############
#load data from folder Brain tissues.
##############
brain_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Brain_tau/Ambrain_log1_sbun_tau_match_sorted_fi.txt", header = TRUE)
str(brain_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_brain.pdf", width=8, height=8)
ggplot(brain_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

####
cor.test(sqrt(abs(brain_tau$logFC.XYbrain.XXrain)),sqrt(brain_tau$tau),
         method = "spearman",
         conf.level = 0.95)

###
data:  sqrt(abs(brain_tau$logFC.XYbrain.XXrain)) and sqrt(brain_tau$tau)
S = 2.6339e+12, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.1517634 
###

###some stats
wilcox.test(brain_tau$tau[brain_tau$bias=='female'],brain_tau$tau[brain_tau$bias=='unbias'],exact = FALSE) 
#W = 1551500, p-value < 2.2e-16

wilcox.test(brain_tau$tau[brain_tau$bias=='male'],brain_tau$tau[brain_tau$bias=='unbias'],exact = FALSE) 
#W = 1235700, p-value = 0.0001819

wilcox.test(brain_tau$tau[brain_tau$bias=='male'],brain_tau$tau[brain_tau$bias=='female'],exact = FALSE) 
#W = 1329, p-value = 1.093e-07

###dNdS and tau in Brain tissues
brain_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Brain_tau/Ambrain_log1_sb_tau_match_dnds_sorted_fi_shareunbias.txt", header = TRUE)
str(brain_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/brain_sbun_dnds_shareunbias.pdf", width=8, height=8)
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
#W = 5540, p-value = 0.1364
wilcox.test(brain_tau_dnds$dNdS[brain_tau_dnds$bias=='male'],brain_tau_dnds$dNdS[brain_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 23444, p-value = 0.1086
wilcox.test(brain_tau_dnds$dNdS[brain_tau_dnds$bias=='female'],brain_tau_dnds$dNdS[brain_tau_dnds$bias=='male'],exact = FALSE) 
#W = 74, p-value = 0.04997

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_brain_shareunbias.pdf", width=8, height=8)
ggplot(brain_tau_dnds, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

###some stats
wilcox.test(brain_tau_dnds$tau[brain_tau_dnds$bias=='female'],brain_tau_dnds$tau[brain_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 11708, p-value = 0.003417
wilcox.test(brain_tau_dnds$tau[brain_tau_dnds$bias=='male'],brain_tau_dnds$tau[brain_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 30195, p-value = 4.479e-06
wilcox.test(brain_tau_dnds$tau[brain_tau_dnds$bias=='female'],brain_tau_dnds$tau[brain_tau_dnds$bias=='male'],exact = FALSE) 
#W = 129, p-value = 0.9859

#with color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_dnds_brain_colors.pdf", width=8, height=8)
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
Df Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)         1  0.465 0.46481 50.9438 1.049e-12 ***
  bias              2  0.036 0.01805  1.9782    0.1384    
sqrt(tau):bias    2  0.022 0.01086  1.1902    0.3042    
Residuals      6829 62.308 0.00912   
##########

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_brain_exp.pdf", width=8, height=8)
p5 <-ggplot2.scatterplot(data=brain_tau, xName='asblogFC.XYbrain.Xxrain',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="Brain: Log2(male/female)", y="Tissue specificity") +
  scale_fill_manual(values = c("grey40"))
dev.off()


##############
#load data from folder Liver tissues.
##############
liver_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Liver_tau/Amliver_log1_sbun_tau_match_sorted.txt", header = TRUE)
str(liver_tau)

abs <- abs(liver_tau$logFC.XYliver.XXliver)
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/tau_amm_liver.pdf", width=8, height=8)
ggplot(liver_tau, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias",labels=c("XX","XY","unbias")) +
  geom_boxplot() +
  ylim(0,1) +
  scale_x_discrete(labels=c("XX", "XY","Unbias"),name="Sex bias") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11))
dev.off()

cor.test(sqrt(abs(liver_tau$logFC.XYliver.XXliver)),sqrt(liver_tau$tau),
         method = "spearman",
         conf.level = 0.95)

###
data:  sqrt(abs(liver_tau$logFC.XYliver.XXliver)) and sqrt(liver_tau$tau)
S = 9.3843e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.1675863 

###

###some stats
wilcox.test(liver_tau$tau[liver_tau$bias=='female'],liver_tau$tau[liver_tau$bias=='unbias'],exact = FALSE) 
#W = 1819400, p-value = 2.897e-16
wilcox.test(liver_tau$tau[liver_tau$bias=='male'],liver_tau$tau[liver_tau$bias=='unbias'],exact = FALSE) 
#W = 1266700, p-value = 7.196e-10
wilcox.test(liver_tau$tau[liver_tau$bias=='female'],liver_tau$tau[liver_tau$bias=='male'],exact = FALSE) 
#W = 7352, p-value = 0.4464

####
liver_tau_dnds <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/Liver_tau/Amliver_log1_sb_shareunbias_tau_dnds_fi.txt", header = TRUE)
str(liver_tau_dnds)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/liver_sbun_dnds.pdf", width=8, height=8)
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
#W = 23308, p-value = 0.05493
wilcox.test(liver_tau_dnds$dNdS[liver_tau_dnds$bias=='male'],liver_tau_dnds$dNdS[liver_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 35469, p-value = 0.002628
wilcox.test(liver_tau_dnds$dNdS[liver_tau_dnds$bias=='female'],liver_tau_dnds$dNdS[liver_tau_dnds$bias=='male'],exact = FALSE) 
#W = 426.5, p-value = 0.7359


pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/liver_sbun_tau_sharebias.pdf", width=8, height=8)
ggplot(liver_tau_dnds, aes(x=bias, y=tau, fill=bias)) + scale_fill_manual(values = c("firebrick2","dodgerblue2","grey40"), name="Sex bias") +
  geom_boxplot() +
  ylim(0,1) +
  labs(x='Sex bias', y='dN/dS') +
  scale_x_discrete(labels=c("XX", "XY", "unbias")) +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###some stats
wilcox.test(liver_tau_dnds$tau[liver_tau_dnds$bias=='female'],liver_tau_dnds$tau[liver_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 29292, p-value = 3.9e-06
wilcox.test(liver_tau_dnds$tau[liver_tau_dnds$bias=='male'],liver_tau_dnds$tau[liver_tau_dnds$bias=='unbias'],exact = FALSE) 
#W = 41192, p-value = 2.598e-07
wilcox.test(liver_tau_dnds$tau[liver_tau_dnds$bias=='female'],liver_tau_dnds$tau[liver_tau_dnds$bias=='male'],exact = FALSE) 
#W = 479, p-value = 0.676

####
qqnorm(sqrt(abs(liver_tau$tau)))
hist(sqrt(abs(liver_tau$tau)))
qqnorm(sqrt(abs(liver_tau_dnds$dNdS)))
hist(sqrt(abs(liver_tau_dnds$dNdS)))

y3 <- lm(sqrt(dNdS) ~ sqrt(tau) * bias, liver_tau_dnds)
anova(y3)

##
Response: sqrt(dNdS)
Df  Sum Sq Mean Sq F value    Pr(>F)    
sqrt(tau)         1  0.4148 0.41482 45.3508 2.299e-11 ***
  bias              2  0.0813 0.04066  4.4456   0.01188 *  
  sqrt(tau):bias    2  0.0564 0.02819  3.0817   0.04616 *  
  Residuals      1579 14.4431 0.00915                     
---
##
#with color
  pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_dnds_liver_colors.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="dNdS", y="Tau", color ="bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_exp_liver_colors.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau, xName=abs,yName='tau', ylim=c(0,1),size=2, groupName='bias',groupColors=c("firebrick4","dodgerblue4","grey50"), addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="abs(logFC.XYliver.XXliver)", y="Tau", color ="bias")
dev.off()

#without color
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_liver_dnds.pdf", width=8, height=8)
ggplot2.scatterplot(data=liver_tau_dnds, xName='dNdS',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="dNdS", y="Tau") +
  scale_fill_manual(values = c("grey40"))
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/tau/scatter_abs_tau_exp_liver_nocolors.pdf", width=8, height=8)
p4 <- ggplot2.scatterplot(data=liver_tau, xName='abslogFC.XYliver.Xxliver',yName='tau', ylim=c(0,1),size=2,addRegLine=TRUE, addConfidenceInterval=TRUE,color='grey40')  +
  labs(x="Liver: Log2(male/female)", y="Tissue specificity") +
  scale_fill_manual(values = c("grey40"))
dev.off()


###
alltissues_tau <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/am_tau_8tissues_fi.txt", header = TRUE)
str(alltissues_tau)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/alltissue_tau.pdf", width=8, height=8)

ggplot(data = alltissues_tau, aes(x=stage, y= tau,fill=stage)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1.1)) +
  scale_fill_manual(values = c("grey","grey","grey","grey","grey","grey","grey","grey")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

fourtissues_tau_top30 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_g2346gonadliver.txt", header = TRUE)
str(fourtissues_tau_top30)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/4tissue_tau_top30.pdf", width=8, height=8)

ggplot(data = fourtissues_tau_top30, aes(x=stage, y= tau,fill=stage)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1.1)) +
  scale_fill_manual(values = c("grey","grey","grey","grey")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

alltissues_tau_top70 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_alltissues_top70.txt", header = TRUE)
str(alltissues_tau_top70)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/alltissue_tau_top70.pdf", width=8, height=8)

ggplot(data = alltissues_tau_top70, aes(x=stage, y= tau,fill=stage)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey","grey","grey","grey","grey","grey","grey","grey")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#g23

g23_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg23.txt", header = TRUE)
str(g23_all)

g23_all_sort <- data.frame(g23_all[order(g23_all$logCPM, decreasing=T),])

nr_g23 <- nrow(g23_all_sort)*0.9
nr_g23

nr_g23 <- nrow(g23_all_sort)*0.7
nr_g23

nr_g23_top30 <- nrow(g23_all_sort)*0.3
nr_g23_top30

nr_g23_top10 <- nrow(g23_all_sort)*0.1
nr_g23_top10

g23_top10 <- data.frame(g23_all_sort[1:2387,])
write.table(g23_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg23_top10.txt", sep="\t", col.names=T)

g23_top90 <- data.frame(g23_all_sort[1:21490,])
write.table(g23_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg23_top90.txt", sep="\t", col.names=T)

g23_top30 <- data.frame(g23_all_sort[1:7163,])
write.table(g23_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg23_top30.txt", sep="\t", col.names=T)

g23_top70 <- data.frame(g23_all_sort[1:16714,])
write.table(g23_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg23_top70.txt", sep="\t", col.names=T)

g23_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amg23_topall.txt", header = TRUE)
str(g23_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_amg23_topall", width=8, height=8)

ggplot(data = g23_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()
##g27

g27_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg27.txt", header = TRUE)
str(g27_all)

g27_all_sort <- data.frame(g27_all[order(g27_all$logCPM, decreasing=T),])

nr_g27 <- nrow(g27_all_sort)*0.9
nr_g27

nr_g27 <- nrow(g27_all_sort)*0.7
nr_g27

nr_g27_top30 <- nrow(g27_all_sort)*0.3
nr_g27_top30

nr_g27 <- nrow(g27_all_sort)*0.1
nr_g27

g27_top10 <- data.frame(g27_all_sort[1:2437,])
write.table(g27_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg27_top10.txt", sep="\t", col.names=T)

g27_top90 <- data.frame(g27_all_sort[1:21490,])
write.table(g27_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg27_top90.txt", sep="\t", col.names=T)

g27_top30 <- data.frame(g27_all_sort[1:7163,])
write.table(g27_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg27_top30.txt", sep="\t", col.names=T)

g27_top70 <- data.frame(g27_all_sort[1:17059,])
write.table(g27_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg27_top70.txt", sep="\t", col.names=T)

g27_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amg27_topall.txt", header = TRUE)
str(g27_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_amg27_topall", width=8, height=8)

ggplot(data = g27_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#g31


g31_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg31.txt", header = TRUE)
str(g31_all)

g31_all_sort <- data.frame(g31_all[order(g31_all$logCPM, decreasing=T),])

nr_g31 <- nrow(g31_all_sort)*0.9
nr_g31

nr_g31 <- nrow(g31_all_sort)*0.7
nr_g31

nr_g31_top30 <- nrow(g31_all_sort)*0.3
nr_g31_top30

nr_g31_top10 <- nrow(g31_all_sort)*0.1
nr_g31_top10

g31_top10 <- data.frame(g31_all_sort[1:2375,])
write.table(g31_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg31_top10.txt", sep="\t", col.names=T)

g31_top90 <- data.frame(g31_all_sort[1:21490,])
write.table(g31_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg31_top90.txt", sep="\t", col.names=T)

g31_top30 <- data.frame(g31_all_sort[1:7163,])
write.table(g31_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg31_top30.txt", sep="\t", col.names=T)

g31_top70 <- data.frame(g31_all_sort[1:16629,])
write.table(g31_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg31_top70.txt", sep="\t", col.names=T)

g31_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amg31_topall.txt", header = TRUE)
str(g31_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_g31_topall", width=8, height=8)

ggplot(data = g31_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#g43
#g43

g43_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg43.txt", header = TRUE)
str(g43_all)

g43_all_sort <- data.frame(g43_all[order(g43_all$logCPM, decreasing=T),])

nr_g43 <- nrow(g43_all_sort)*0.9
nr_g43

nr_g43 <- nrow(g43_all_sort)*0.7
nr_g43

nr_g43_top30 <- nrow(g43_all_sort)*0.3
nr_g43_top30

nr_g43_top10 <- nrow(g43_all_sort)*0.1
nr_g43_top10

g43_top10 <- data.frame(g43_all_sort[1:2557,])
write.table(g43_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg43_top10.txt", sep="\t", col.names=T)

g43_top90 <- data.frame(g43_all_sort[1:21490,])
write.table(g43_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg43_top90.txt", sep="\t", col.names=T)

g43_top30 <- data.frame(g43_all_sort[1:7163,])
write.table(g43_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg43_top30.txt", sep="\t", col.names=T)

g43_top70 <- data.frame(g43_all_sort[1:17091,])
write.table(g43_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg43_top70.txt", sep="\t", col.names=T)

g43_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amg43_topall.txt", header = TRUE)
str(g43_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_g43_topall", width=8, height=8)

ggplot(data = g43_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#brain
#brain

brain_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Ambrain.txt", header = TRUE)
str(brain_all)

brain_all_sort <- data.frame(brain_all[order(brain_all$logCPM, decreasing=T),])

nr_brain <- nrow(brain_all_sort)*0.9
nr_brain

nr_brain <- nrow(brain_all_sort)*0.7
nr_brain

nr_brain_top30 <- nrow(brain_all_sort)*0.3
nr_brain_top30

nr_brain_top10 <- nrow(brain_all_sort)*0.1
nr_brain_top10

brain_top10 <- data.frame(brain_all_sort[1:2655,])
write.table(brain_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/ambrain_top10.txt", sep="\t", col.names=T)

brain_top90 <- data.frame(brain_all_sort[1:23899,])
write.table(brain_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/ambrain_top90.txt", sep="\t", col.names=T)

brain_top30 <- data.frame(brain_all_sort[1:7163,])
write.table(brain_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/ambrain_top30.txt", sep="\t", col.names=T)

brain_top70 <- data.frame(brain_all_sort[1:18588,])
write.table(brain_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/ambrain_top70.txt", sep="\t", col.names=T)

brain_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_ambrain_topall.txt", header = TRUE)
str(brain_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_brain_topall", width=8, height=8)

ggplot(data = brain_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()
##g46

g46_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg46.txt", header = TRUE)
str(g46_all)

g46_all_sort <- data.frame(g46_all[order(g46_all$logCPM, decreasing=T),])

nr_g46 <- nrow(g46_all_sort)*0.9
nr_g46

nr_g46 <- nrow(g46_all_sort)*0.7
nr_g46

nr_g46_top30 <- nrow(g46_all_sort)*0.3
nr_g46_top30
nr_g46_top10 <- nrow(g46_all_sort)*0.1
nr_g46_top10

g46_top10 <- data.frame(g46_all_sort[1:2494,])
write.table(g46_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg46_top10.txt", sep="\t", col.names=T)

g46_top90 <- data.frame(g46_all_sort[1:21490,])
write.table(g46_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg46_top90.txt", sep="\t", col.names=T)

g46_top30 <- data.frame(g46_all_sort[1:7482,])
write.table(g46_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg46_top30.txt", sep="\t", col.names=T)

g46_top70 <- data.frame(g46_all_sort[1:17458,])
write.table(g46_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amg46_top70.txt", sep="\t", col.names=T)


g46_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amg46_topall.txt", header = TRUE)
str(g46_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_g46_topall", width=8, height=8)

ggplot(data = g46_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#gonad

gonad_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amgonad.txt", header = TRUE)
str(gonad_all)

gonad_all_sort <- data.frame(gonad_all[order(gonad_all$logCPM, decreasing=T),])

nr_gonad <- nrow(gonad_all_sort)*0.9
nr_gonad

nr_gonad_top30 <- nrow(gonad_all_sort)*0.3
nr_gonad_top30

nr_gonad_top70 <- nrow(gonad_all_sort)*0.7
nr_gonad_top70
nr_gonad_top10 <- nrow(gonad_all_sort)*0.1
nr_gonad_top10

gonad_top10 <- data.frame(gonad_all_sort[1:2256,])
write.table(gonad_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amgonad_top10.txt", sep="\t", col.names=T)

gonad_top90 <- data.frame(gonad_all_sort[1:21490,])
write.table(gonad_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amgonad_top90.txt", sep="\t", col.names=T)

gonad_top30 <- data.frame(gonad_all_sort[1:6768,])
write.table(gonad_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amgonad_top30.txt", sep="\t", col.names=T)

gonad_top70 <- data.frame(gonad_all_sort[1:15792,])
write.table(gonad_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amgonad_top70.txt", sep="\t", col.names=T)

gonad_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amgonad_topall.txt", header = TRUE)
str(gonad_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_gonad_topall", width=8, height=8)

ggplot(data = gonad_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()

#liver
##

liver_all<- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amliver.txt", header = TRUE)
str(liver_all)

liver_all_sort <- data.frame(liver_all[order(liver_all$logCPM, decreasing=T),])

nr_liver <- nrow(liver_all_sort)*0.9
nr_liver

nr_liver <- nrow(liver_all_sort)*0.7
nr_liver

nr_liver_top30 <- nrow(liver_all_sort)*0.3
nr_liver_top30

nr_liver_top10 <- nrow(liver_all_sort)*0.1
nr_liver_top10

liver_top10 <- data.frame(liver_all_sort[1:1892,])
write.table(liver_top10, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amliver_top10.txt", sep="\t", col.names=T)

liver_top90 <- data.frame(liver_all_sort[1:21490,])
write.table(liver_top90, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amliver_top90.txt", sep="\t", col.names=T)

liver_top30 <- data.frame(liver_all_sort[1:5677,])
write.table(liver_top30, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amliver_top30.txt", sep="\t", col.names=T)
liver_top70 <- data.frame(liver_all_sort[1:13247,])
write.table(liver_top70, "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/amliver_top70.txt", sep="\t", col.names=T)

liver_tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amliver_topall.txt", header = TRUE)
str(liver_tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_liver_topall", width=8, height=8)

ggplot(data = liver_tau_topall, aes(x=perc, y= tau,fill=perc)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,width=0.7,alpha=0.9) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("grey80","grey60","grey40","grey20")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.position="none")
dev.off()


###for all 8 tissues.

tau_topall <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/tau_amall_topall_fifi.txt", header = TRUE)
str(tau_topall)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/tau_alltissues_topall.pdf", width=8, height=8)

ggplot(data = tau_topall, aes(x=perc, y= tau,fill=stage)) + 
  geom_boxplot(notch=TRUE,outlier.shape=NA,position=position_dodge(0.7), width=0.55,alpha=0.8) +
  coord_cartesian(ylim = c(0,0.3,0.6,0.9,1)) +
  scale_fill_manual(values = c("steelblue1","steelblue3","steelblue4","purple1","purple4","springgreen1","springgreen3","springgreen4")) +
  theme_set(theme_bw(base_size=12)) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) 

dev.off()

# G23 unbias correlation with tau

tau_unbias_amg23 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg23_tau.txt", header = TRUE)
str(tau_unbias_amg23)

qqnorm(sqrt(tau_unbias_amg23$tau))
hist(sqrt(tau_unbias_amg23$tau))

cor.test(sqrt(tau_unbias_amg23$tau),sqrt(tau_unbias_amg23$logFC.XY23.XX23),
         method = "spearman",
         conf.level = 0.95)

###
S = 2.029e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.225896 
###

tau_unbias_amg27 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg27_tau.txt", header = TRUE)
str(tau_unbias_amg27)

qqnorm(sqrt(tau_unbias_amg27$tau))
hist(sqrt(tau_unbias_amg27$tau))

cor.test(sqrt(tau_unbias_amg27$tau),sqrt(tau_unbias_amg27$logFC.XY27.XX27),
         method = "spearman",
         conf.level = 0.95)
###
S = 2.312e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2233488 
###
tau_unbias_amg31 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/sexsp_tau/cpm_tau/LogCPM_0.05_Amg31_tau.txt", header = TRUE)
str(tau_unbias_amg31)

qqnorm(sqrt(tau_unbias_amg31$tau))
hist(sqrt(tau_unbias_amg31$tau))

cor.test(sqrt(tau_unbias_amg31$tau),sqrt(tau_unbias_amg31$logFC.XY31.XX31),
         method = "spearman",
         conf.level = 0.95)
###
S = 2.399e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.2247824 
###

pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/correlation_5tissuesv2.pdf")

ggarrange(p1, p2, p5, p3, p4, labels = c("A","B","C","D","E"),
          ncol = 2, nrow = 3)

dev.off()
