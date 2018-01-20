#intall R packages and loading libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)


#set up the working directory
setwd("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds")
#results directoy
setwd("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

dn_ds_all<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds/all_dnds_exp_sorted.txt", header = T)
str(dn_ds_all)

sex_auto<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds/sexchr_auto_dnds_all.txt", header = T)
str(sex_auto)
head(sex_auto)

#dn/ds plot for Faster-X
###all chromosome separately
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_fasterX_10chr.pdf", width=10)
ggplot(sex_auto, aes(x=chrid, y=dNdS, fill=chrid)) + 
  scale_fill_manual(values = c("firebrick2","firebrick2","grey","grey","grey","grey","grey","grey","grey","grey")) +
  theme(legend.position="none") +
  geom_boxplot() +
  ylim(0,0.5) +
  labs(x='Chromosome', y='dn/ds') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

#sex chromosome vs autosomes
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_fasterX.pdf")
ggplot(sex_auto, aes(x=chrtype, y=dNdS, fill=chrtype)) + 
  scale_fill_manual(values = c("grey","firebrick2")) +
  theme(legend.position="none") +
  geom_boxplot() +
  ylim(0,0.5) +
  labs(x='Chromosome', y='dn/ds') +
  scale_x_discrete(labels = c("Autosome","Sex chromosome")) +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

###if separate to compare Chr01, Chr02 and autosome
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_fasterX_chr0102.pdf")
ggplot(sex_auto, aes(x=twosexchr, y=dNdS, fill=twosexchr)) + 
  scale_fill_manual(values = c("grey","firebrick2","firebrick2")) +
  theme(legend.position="none") +
  geom_boxplot() +
  ylim(0,0.5) +
  labs(x='Chromosome', y='dn/ds') +
  scale_x_discrete(labels = c("Autosome","Chr01","Chr02")) +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

mean(sex_auto$dNdS[sex_auto$chrid=='Chr01'])  #0.08367571
mean(sex_auto$dNdS[sex_auto$chrid=='Chr02']) # 0.08482178
mean(sex_auto$dNdS[sex_auto$chrtype=='Autosome'])  #0.08202889

wilcox.test(sex_auto$dNdS[sex_auto$twosexchr=='Autosome'], sex_auto$dNdS[sex_auto$twosexchr=='Chr01'])
#W = 2745200, p-value = 0.2957
wilcox.test(sex_auto$dNdS[sex_auto$twosexchr=='Chr02'], sex_auto$dNdS[sex_auto$twosexchr=='Autosome'])
#W = 2453400, p-value = 0.06201
wilcox.test(sex_auto$dNdS[sex_auto$twosexchr=='Chr01'], sex_auto$dNdS[sex_auto$twosexchr=='Chr02'])
#data:  sex_auto$dNdS[sex_auto$chr == "Chr01"] and sex_auto$dNdS[sex_auto$chr == "Chr02"]
#W = 566390, p-value = 0.4656

#boxplot
#dn/ds
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_dnds_all.pdf", width=10)
ggplot(dn_ds_all, aes(x=sexbias, y=dnds,fill=(sexbias))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("F", "M", "U"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,0.6)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#ds
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_ds_all.pdf", width=10)
ggplot(dn_ds_all, aes(x=sexbias, y=ds,fill=(sexbias))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("F", "M", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dS", limits = c(0,2)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#dn
pdf(file="/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/amm_dn_all.pdf", width=10)
ggplot(dn_ds_all, aes(x=sexbias, y=dn,fill=(sexbias))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("F", "M", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN", limits = c(0,1)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

G31 <- subset(dn_ds_all, dn_ds_all$stage=='G31')
G43 <- subset(dn_ds_all, dn_ds_all$stage=='G43')
G46 <- subset(dn_ds_all, dn_ds_all$stage=='G46')
gonad <- subset(dn_ds_all, dn_ds_all$stage=='Gonad')
liver <- subset(dn_ds_all, dn_ds_all$stage=='Liver')
brain <- subset(dn_ds_all, dn_ds_all$stage=='Brain')


#G43
wilcox.test(G43$dnds[G43$sexbias=='female'], G43$dnds[G43$sexbias=='male']) # W = 60, p-value = 0.9149
wilcox.test(G43$dnds[G43$sexbias=='female'], G43$dnds[G43$sexbias=='unbias']) #W = 202960, p-value = 6.255e-06
wilcox.test(G43$dnds[G43$sexbias=='male'], G43$dnds[G43$sexbias=='unbias']) #W = 14262, p-value = 0.2533. This is not correct, as only 3 male values.

#G46
wilcox.test(G46$dnds[G46$sexbias=='male'], G46$dnds[G46$sexbias=='unbias']) #W = 441010, p-value = 8.044e-08
wilcox.test(G46$dnds[G46$sexbias=='female'], G46$dnds[G46$sexbias=='unbias']) #W = 4095500, p-value = 0.1994
wilcox.test(G46$dnds[G46$sexbias=='female'], G46$dnds[G46$sexbias=='male']) #W = 95099, p-value = 3.095e-06

#gonad
wilcox.test(gonad$dnds[gonad$sexbias=='male'], gonad$dnds[gonad$sexbias=='unbias']) #W = 1855500, p-value = 0.06863
wilcox.test(gonad$dnds[gonad$sexbias=='female'], gonad$dnds[gonad$sexbias=='unbias']) #W = 2430100, p-value = 0.002517
wilcox.test(gonad$dnds[gonad$sexbias=='female'], gonad$dnds[gonad$sexbias=='male']) #W = 1349300, p-value = 0.3569

#brain
wilcox.test(brain$dnds[brain$sexbias=='male'], brain$dnds[brain$sexbias=='unbias']) #W = 103310, p-value = 0.1372
wilcox.test(brain$dnds[brain$sexbias=='female'], brain$dnds[brain$sexbias=='unbias']) #W = 24332, p-value = 0.1197
wilcox.test(brain$dnds[brain$sexbias=='female'], brain$dnds[brain$sexbias=='male']) #W = 74, p-value = 0.04881

#liver
wilcox.test(liver$dnds[liver$sexbias=='male'],liver$dnds[liver$sexbias=='unbias']) #W = 123380, p-value = 0.005197
wilcox.test(liver$dnds[liver$sexbias=='female'],liver$dnds[liver$sexbias=='unbias']) #W = 81166, p-value = 0.0796
wilcox.test(liver$dnds[liver$sexbias=='female'],liver$dnds[liver$sexbias=='male']) #W = 426.5, p-value = 0.7359

#combined all sex biased in three stages.
ggplot2.boxplot(data=dn_ds_all, xName='sb', yName='dnds', group='sb', groupColors=c("firebrick2","dodgerblue2","grey"), ylim=c(0,0.1))
ggplot2.boxplot(data=dn_ds_all, xName='sb', yName='dnds', group='sb', groupColors=c("firebrick2","dodgerblue2","grey"), ylim=c(0,1.1))
ggplot(share_dnds_4346, aes(x=sb, y=dnds,fill=sb)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  scale_x_discrete(labels=c("XX", "XY'", "Unbias"),name="Sex bias") +
  labs(y='dN/dS') +
  theme(legend.position="none") + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='female'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) # W=180, P=0.29
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='female'], share_dnds_4346$dnds[share_dnds_4346$sb=='unbias']) # W=46422, P=2.7e-05
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='male'], share_dnds_4346$dnds[share_dnds_4346$sb=='unbias']) # W=11690, P=0.0016

wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='M-bias']) # W857500, P=0.40
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='unbias']) # W85223700, P=0.50
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='M-bias'], dn_ds_all$dnds[dn_ds_all$sb=='unbias']) # W5255400, P=0.64

p <- ggplot(dn_ds_4346, aes(x=sb, y=dnds,color=sb)) + 
  geom_boxplot(notch = FALSE)
p + labs(x="Sex bias", y="dN/dS", aex.lab=2, aex.axis=1.5, aex.main=2) +
  theme(legend.position="none")
wilcox.test(dn_ds_4346$dnds[dn_ds_all$sb=='unbias'], dn_ds_all$dnds[dn_ds_all$sb=='F-bias'])

# G43 46 stages, all categories of sex bias and unbiased genes
ggplot(dn_ds_4346_allcat, aes(x=sb, y=dnds, fill=interaction(factor(stage),sb))) + geom_boxplot(notch = FALSE) +
  theme(legend.position="none") +
scale_fill_manual(values = c("firebrick2","dodgerblue2","grey90","firebrick4","dodgerblue4","grey50")) +
  #geom_rect(xmin = 1.5, xmax = 10.5, ymin=-Inf, ymax=Inf,alpha=0.5, fill="grey75") +
  geom_vline(aes(xintercept=1.5), linetype="blank") +
  geom_vline(aes(xintercept=10.5), linetype="blank") +
  labs(x='Sex bias', y='dn/ds')

wilcox.test(dn_ds_4346$dnds[dn_ds_all$sb=='unbias'], dn_ds_all$dnds[dn_ds_all$sb=='F-bias'])

#G46 log2>=2 

ggplot(dnds_46_log2, aes(x=sb, y=dnds, fill=sb)) + geom_boxplot(notch = FALSE) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey90")) +
  labs(x='Sex bias at G46 (Log2 >=2)', y='dn/ds')

wilcox.test(dnds_46_log2$dnds[dnds_46_log2$sb=='unbias'], dnds_46_log2$dnds[dnds_46_log2$sb=='F-bias'])

#G46 Log2>=3
ggplot(dnds_46_log3, aes(x=sb, y=dnds, fill=sb)) + geom_boxplot(notch = FALSE) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey90")) +
  labs(x='Sex bias at G46 (Log2 >=2)', y='dn/ds')

wilcox.test(dnds_46_log3$dnds[dnds_46_log3$sb=='F-bias'], dnds_46_log3$dnds[dnds_46_log3$sb=='M-bias'])

#scatter plot
ggplot2.scatterplot(data=share_dnds_abs4346, xName='X43_XY_XX_abs',yName='dnds', 
                    #groupName='sb', size=4,
                    #groupColors=c("firebrick4","dodgerblue4","grey50"),
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) gene expression at G43", y="dn/ds", color ="Sex bias") +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

ggplot2.scatterplot(data=share_dnds_abs4346, xName='X46_XY_XX_abs',yName='dnds', 
                    addRegLine=TRUE, addConfidenceInterval=TRUE) +
                    #groupName='sb', size=4,
                    #groupColors=c("firebrick4","dodgerblue4","grey50")) +
  labs(x="Log2 (XY/XX) gene expression at G46", y="dn/ds", color ="Sex bias") +
    theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
    theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

ggplot(data=share_dnds_4346, aes(x=sb, y=dnds, fill=sb)) + geom_boxplot(notch = FALSE) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey90")) +
  labs(x='Shared expressed genes at G43 and G46', y='dn/ds') +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=11),axis.text.y = element_text(colour="black",size=11)) +
  scale_x_discrete(labels = c("Female", "Male", "Unbias"))

y<- glm(dnds~X43_XY_XX_abs,data=share_dnds_abs4346,family=quasibinomial)
anova(y)
summary(y)

y1<- glm(dnds~X46_XY_XX_abs,data=share_dnds_abs4346,family=quasibinomial)
anova(y1)
summary(y1)
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$ID=='46_unbias'], share_dnds_4346$dnds[share_dnds_4346$ID=='share_Fbias'])
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$ID=='46_unbias'], share_dnds_4346$dnds[share_dnds_4346$ID=='share_Mbias'])
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$ID=='share_Mbias'], share_dnds_4346$dnds[share_dnds_4346$ID=='share_Fbias'])

ggplot2.scatterplot(data=share_dnds_4346, xName='X46_XY_XX',yName='dnds', 
                    groupName='sb', size=4,
                    groupColors=c("firebrick4","dodgerblue4","grey50"),
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")

ggplot2.scatterplot(data=share_dnds_abs4346, xName='X46_XY_XX_abs',yName='dnds', 
                  addRegLine=TRUE, addConfidenceInterval=TRUE,
                  groupName='sb', size=4,
groupColors=c("firebrick4","dodgerblue4","grey50")) +
labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")

y <-lm(data=share_dnds_4346,dnds~X43_XY_XX*sb)
anova(y)
y1 <-lm(data=share_dnds_4346,dnds~X46_XY_XX*sb)
anova(y1)
  
#data is unlikely to be normal istribution, hence it is not proper to use lm test.
data_43 <- data.frame(dn_ds_all$sb[dn_ds_all$stage=='43'],dn_ds_all$dnds[dn_ds_all$stage=='43'])
lm_43 <- lm(dn_ds_all$dnds[dn_ds_all$stage=='43']~dn_ds_all$sb[dn_ds_all$stage=='43'], data_43)
anova(lm_43)
summary(lm_43)
#G43, unbias vs. F-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='43_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='43_F-bias'])

########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "43_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_F-bias"]
W = 65316, p-value = 2.357e-06
alternative hypothesis: true location shift is not equal to 0
#########

#G43, unbias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='43_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='43_M-bias'])

#########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "43_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_M-bias"]
W = 11844, p-value = 0.0001248
alternative hypothesis: true location shift is not equal to 0
#########

#G43, F-bias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='43_F-bias'], dn_ds_all$dnds[dn_ds_all$ID=='43_M-bias'])

########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "43_F-bias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_M-bias"]
W = 280, p-value = 0.2222
alternative hypothesis: true location shift is not equal to 0
##############

#G46, F-bias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_F-bias'], dn_ds_all$dnds[dn_ds_all$ID=='46_M-bias'])

########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "46_F-bias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_M-bias"]
W = 816360, p-value = 0.3217
alternative hypothesis: true location shift is not equal to 0
########

#G46, unbias vs. F-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='46_F-bias'])

#########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "46_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_F-bias"]
W = 1327900, p-value = 0.03341
alternative hypothesis: true location shift is not equal to 0
#########

#G46, unbias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='46_M-bias'])

#########
Wilcoxon rank sum test with continuity correction

data:  dn_ds_all$dnds[dn_ds_all$ID == "46_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_M-bias"]
W = 828440, p-value = 0.3571
alternative hypothesis: true location shift is not equal to 0
#########

p1 <- ggplot(dn_ds_all, aes(x=ID, y=dn_ds_all$dn,color= dn_ds_all$sb)) + 
  geom_boxplot(notch = FALSE)
p1 + labs(title="dN across developmental stages",x="Developmental stage", y="dN", aex.lab=2, aex.axis=1.5, aex.main=2)

p2 <- ggplot(dn_ds_all, aes(x=ID, y=dn_ds_all$ds,color= dn_ds_all$sb)) + 
  geom_boxplot(notch = FALSE)
p2 + labs(title="dS across developmental stages",x="Developmental stage", y="dS", aex.lab=2, aex.axis=1.5, aex.main=2)

m0<- glm(dn_ds_all$dnds~dn_ds_all$ID)

violinData <- data.frame(subset(restab_logCPM, grepl("Chr", chr))["chr"], subset(restab_logCPM, grepl("Chr", chr))[paste('logFC.',colnames(cmat)[k], sep="")])
p = ggplot(violinData, aes_string(factor(violinData$chr),paste('logFC.',colnames(cmat)[k], sep="")), ylim)
p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC.',colnames(cmat)[k], sep=""), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
ggsave(file.path(outpath, paste('violin_per_chr_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)


restab_logCPM_subset <- subset(restab_logCPM, restab_logCPM[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use) 
violinData <- data.frame(subset(restab_logCPM_subset, grepl("Chr", chr))["chr"], subset(restab_logCPM_subset, grepl("Chr", chr))[paste('logFC.',colnames(cmat)[k], sep="")])
p = ggplot(violinData, aes_string(factor(violinData$chr),paste('logFC.',colnames(cmat)[k], sep="")), ylim)
p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC.',colnames(cmat)[k],' DE', sep=""), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
ggsave(file.path(outpath, paste('violin_per_chr_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)

}

summary_all <- ddply(dn_ds_all,c("stage","sb"),summarise,mean=mean(dnds),N=length(dnds),sd =sd(dnds),se=sd/sqrt(N))          
summary_all

#ggplot
pd <- position_dodge(0.1)
ggplot(summary_all,aes(x=stage,y=mean,colour=sb, group=sb)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5,position=pd) +
  geom_point(position=pd,size=3)


#barplot
summary_all2 <- summary_all
summary_all2$stage <-factor(summary_all$stage) 
ggplot(summary_all2,aes(x=stage,y=mean,fill=sb)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,position=position_dodge(.9)) +
  xlab("Developmental stage (Gosner)") +
  ylab("Divergence (dN/dS)") +
  scale_fill_hue(name="Sex biased",breaks=c("down","un","up"),labels=c("female-biased","unbiased","male-biased")) +
  ggtitle("Divergence of sex biased and unbiased genes") +
  ylim(0,0.3)
  #scale_y_continuous(breaks = 0:20*4) +
  theme_bw()
  
library(ggplot2)  
sex_bias <-read.table("/Users/Wen-Juan/git/rana_transcriptome/input/sexbias.txt", header = T)
str(sex_bias)
sex_bias <- data.frame(sex_bias)
plot(c(23, 46), c(0, 8500), axes=F, lwd=2, xlab="Developmental stages", ylab="Number", type = "null", col="gray50", main="Sex-biased genes across developmental stages", cex.axis=1.5, cex.lab=1.5, cex.main=2)
axis(1,c(0, 30000000,60000000,90000000, 120000000, 150000000, 180000000, 200000000))
axis(2, c(0, 0.5, 1, 1.5, 2, 2.5))
p <- ggplot(sex_bias, aes(stage,bias)) + geom_bar(aes(fill=sex),
                                             width=0.4, position = position_dodge(width=0.5), stat="identity")
p + labs(title="Sex-biased genes across developmental stages",x="Developmental stage", y="Number", aex.lab=2, aex.axis=1.5, aex.main=2) +axis(1,c(23,27,31,43,46))
+axis(2, c(0, 20, 50, 200, 2000, 5000, 8500))
