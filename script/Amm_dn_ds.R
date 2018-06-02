#intall R packages and loading libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)


#set up the working directory
setwd("~/Rana_Transcriptome/input/dnds")
#results directoy
setwd("~/Rana_Transcriptome/output/figures/")

dn_ds_all_shareunbias<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/5tissues_sb_shareunbias_dnds_sorted.txt", header = T)
str(dn_ds_all_shareunbias)

tau_all<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/5tissues_sbun_log1_tau.txt", header = T)
str(tau_all)

sex_auto_mel<-read.table("/Users/Wen-Juan/my_postdoc/postdoc_manuscripts/Amm_RNAseq/Mel/FasterXY.Results/all_dnds.txt", header = T)
str(sex_auto_mel)


#dn/ds
pdf(file="~/Rana_Transcriptome/output/figures/amm_dnds_all_shareunbias.pdf", width=10)

ggplot(dn_ds_all_shareunbias, aes(x=bias, y=dNdS,fill=(bias))) + 
  geom_boxplot(notch = TRUE,outlier.shape=NA) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("F", "M", "U"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,0.3)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

#dn/ds of sex chromosomes vs. autosomes
pdf(file="~/Rana_Transcriptome/output/figures/amm_dnds_all_meldata.pdf", width=10)
ggplot(sex_auto_mel, aes(x=chrom, y=dNdS,fill=(chrom))) + 
  geom_boxplot(notch = TRUE,outlier.shape=NA) +
  scale_fill_manual(values = c("grey","firebrick2","firebrick2")) +
  theme(legend.position="none") +
  scale_y_continuous(name = "dN/dS", limits = c(0,0.3)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

wilcox.test(sex_auto_mel$dNdS[sex_auto_mel$chrom=='Chr02'], sex_auto_mel$dNdS[sex_auto_mel$chrom=='Auto']) #W = 2052300, p-value = 0.8509
wilcox.test(sex_auto_mel$dNdS[sex_auto_mel$chrom=='Chr01'], sex_auto_mel$dNdS[sex_auto_mel$chrom=='Auto']) #W = 2559000, p-value = 0.8916
wilcox.test(sex_auto_mel$dNdS[sex_auto_mel$chrom=='Chr02'], sex_auto_mel$dNdS[sex_auto_mel$chrom=='Chr01']) #W = 476360, p-value = 0.9635


##tau
pdf(file="~/Rana_Transcriptome/output/figures/amm_tau_all_new.pdf", width=10)

ggplot(tau_all, aes(x=bias, y=tau,fill=(bias))) + 
  geom_boxplot(notch = TRUE,outlier.shape=NA) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("F", "M", "U"),name="Sex bias") +
  scale_y_continuous(name = "tau", limits = c(0,1)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

