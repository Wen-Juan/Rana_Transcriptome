#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#load the corresponding data files.

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg43/'
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amg43/LogCPM_0.05_Amg43 copy.txt",header = T)
str(kdata)

###restrict the analysis to sex-biaxed genes
#kdata <- subset(kdata, kdata$logFC.XY43.XX43<=-1)
#kdata <- subset(kdata, kdata$logFC.XY43.XX43>=1)
#kdata$logFC.XY43.XX43

map.data <- subset(kdata, kdata$start!='NA')
chr1.1.data <- data.frame(subset(map.data, map.data$chr=='Chr01')) 
chr1.2.data <- data.frame(subset(map.data, map.data$chr=='Chr02')) 
chr1.data <- rbind(chr1.1.data,chr1.2.data)
chr2.1.data <- subset(map.data, map.data$chr!='Chr01') 
chr2.data<- subset(chr2.1.data, chr2.1.data$chr!='Chr02') 

#group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY male group at G46
group1.expr.minus <- (map.data$Am1_434 + map.data$Am2_434 + map.data$Am4_434) / 3 #XY male group at G43
#group1.expr.minus <- (map.data$A15MT1 + map.data$A17MT1 + map.data$A6MT1 + map.data$A8MT1 + map.data$A12MT1) / 5 #XY male group in gonad
#group1.expr.minus <- (map.data$A10MB + map.data$A15MB + map.data$A16MB + map.data$A17MB + map.data$A8MB) /5 #XY in brain
#group1.expr.minus <- (map.data$A12ML1 + map.data$A17ML1 + map.data$A15ML1 +map.data$A8ML1)/4 #XY linver

#group2.expr.minus <- (map.data$Am2_464+map.data$Am6_464)/2 #XX male group at G46
group2.expr.minus <- (map.data$Am2_433 + map.data$Am4_435 + map.data$Am5_433) / 3 #XX female group at G43
#group2.expr.minus <- (map.data$A10FO1 + map.data$A17FO1 + map.data$A2FO2 + map.data$A8FO2 + map.data$A6FO1) / 5 #XY male group in gonad
#group2.expr.minus <- (map.data$A10FB + map.data$A12FB + map.data$A15FB + map.data$A16FB + map.data$A17FB) /5 #XX in brain
#group2.expr.minus <- (map.data$A12FL1 + map.data$A17FL1 + map.data$A2FL1 + map.data$A7FL1 +map.data$A8FL1)/5 #XX linver

group1.expr <- group1.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + min(abs(c(group1.expr.minus, group2.expr.minus)))


map.data$ratio <- log2(group1.expr/group2.expr)
map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

#ggplot of boxplot
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/log2ratio_liver.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(-8,8)) + 
  geom_boxplot() +
  theme(legend.position="none") +
  labs(x='Chromosome', y='Log2 ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))
dev.off()

###G46 no sex reversals

####
##
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 269820, p-value = 0.6143
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1734400, p-value = 0.9119
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1565500, p-value = 0.453
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 3512100, p-value = 0.4876

###
###G46 no sex reversals

####
##G43
##
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 293860, p-value = 0.2822
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1972500, p-value = 0.8631
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1796200, p-value = 0.1252
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 4026500, p-value = 0.2345

#gonad

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 183190, p-value = 0.164
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1113200, p-value = 0.6068
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 951590, p-value = 0.1998
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 2140700, p-value = 0.5619

##brain

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 329790, p-value = 0.8941
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 2161000, p-value = 0.4588
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1877400, p-value = 0.5987
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 4266700, p-value = 0.3534

##liver

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'],exact = FALSE) 
##W = 159060, p-value = 0.8283
wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 1036600, p-value = 0.7458
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
##W = 853180, p-value = 0.5964
t <- rbind(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr=='Chr02'])
wilcox.test(t,map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 
#W = 2013500, p-value = 0.4988

#calculate means by sliding windows
install.packages("RcppRoll")
library("RcppRoll")
install.packages("zoo")
library("zoo")
install.packages("microbenchmark")
library("microbenchmark")
require(zoo)


dn_ds <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/dnds/amm88perc_annotation_dnds.txt", header = T)
str(dn_ds)

map.data <- dn_ds
chr1_1_start <-map.data$start[map.data$chr=='Chr01']
chr_1_end <- max(map.data$start[map.data$chr=='Chr01'])

zoo.dat <- zoo(map.data$dnds[map.data$chr=='Chr01'], c(chr1_1_start,chr_1_end))
y <- rollapply(zoo.dat, 50, FUN = mean, align = 'center', na.rm=TRUE) 

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_9_dnds_sexchr.pdf", width=8, height=8)
par(mar=c(5,5,4,3)+0.6) 
#plot(c(0, chr_1_end), c(0,2), axes=F, lwd=2, xlab="Position (bp)", ylab="gene expression ratio log2(XY0/XX)", cex.axis=1.5, cex.lab=1.2, col="white")
plot(c(0, chr_1_end), c(0,1.5), axes=F, lwd=2, xlab="Position (bp)", ylab="dN/dS", cex.axis=1.5, cex.lab=1.2, col="white")
axis(1,c(0, 30000000,60000000,90000000, 120000000, 150000000, 180000000, 200000000))
axis(2, c(0,0.2,0.4,0.6,0.8,1,1.2,1.5))
points(map.data$start[map.data$chr=='Chr01'], map.data$dnds[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="gray50", main="",cex.axis=1.5)
#points(DC.data1dmrt$start[DC.data1dmrt$chr=='Chr01'],DC.data1dmrt$dN_dS[DC.data1dmrt$chr=='Chr01'], pch=20, lwd=3, type="p",col="red", main="",cex.axis=1.5)
lines(y,col="blue",lwd=4)
abline(v=116440117, col="blue",lwd=2,lty=3)
dev.off()

