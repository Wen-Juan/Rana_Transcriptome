## need manual changing
compname <- 'Amm46male' #'tv46' #G46
contrast.name <- 'logFC.M46.SR46' # logFC.XYm.XXmt
group1 <- 'XY'
group2 <- 'XX'

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amm46male/'
outpath <- ("~/my_postdoc/useful_scripts/Rana_Transcriptome/output/")
kdata <- read.table(file.path(datapath, compname, paste('LogCPM_0.05_', compname, '.txt', sep="")), header=T)
kdata <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amm46male/LogCPM_0.05_Amm46male_v2.txt",header = T)
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv23_nombias.txt", header=T) #G23 all genes on Chr except male biased genes
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv27_nombias.txt", header=T) #G27 all genes on Chr except male biased genes
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv31_nombias.txt", header=T) #G31 all genes on Chr except male biased genes
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv31_noauto.txt", header=T) #G31 all genes on Chr except male biased genes
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv43_nombias.txt", header=T) #G43 all genes on Chr except male biased genes
kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/output/SB/LogCPM_0.05_tv46_nombias.txt", header=T) #G46 all genes on Chr except male biased genes

# kdata <- read.table("/Users/Wen-Juan/git/rana_transcriptome/input/sexchr_transc_finallist_dnds.txt", header=T)
str(kdata)

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- subset(map.data, map.data$chr=='Chr01') 
chr2.data <- subset(map.data, map.data$chr!='Chr01') 

## need chancing, look at design
group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY group #G46
group2.expr.minus <- map.data$Am5_461 #XX group #G46


group1.expr <- group1.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- group1.expr/group2.expr

map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

write.table(map.data, "~/my_postdoc/useful_scripts/Rana_Transcriptome/output/map.data", sep="\t")

par(mfrow=c(2,2)) 
par(mar=c(5,4,4,2))
dev.copy(pdf, file.path(datapath, compname, paste('DC_', contrast.name, '.pdf', sep="")), width=12, height=18)
plot(density(chr2.data[[contrast.name]]), xlim=range( c(chr1.data[[contrast.name]], chr2.data[[contrast.name]]) ), cex.main=1.4, cex.lab=1.2, main='Density plot') #ylim= c(0,2.5),
lines(density(chr1.data[[contrast.name]]), col=2)
#meanchr01 <- mean(chr1.data[[contrast.name]])
#meanchr02 <- mean(chr2.data[[contrast.name]])
abline(v= -0.55, col=1, lwd=1, lty=2)
abline(v= -0.9, col=2, lwd=1, lty=2)
abline(v= 0.85, col=1, lwd=1, lty=2)
abline(v= 0.65, col=2, lwd=1, lty=2)
rug(chr2.data[[contrast.name]], col=1, ticksize=0.01, line=2.5)
rug(chr1.data[[contrast.name]], col=2, ticksize=0.01, line=3.0)
legend('topleft', inset=0.05, legend=c('Autosomes', 'Chr 1'), pch =c(15,15), col=c(1,2), cex=0.8 )

par(mar=c(5,5,4,4))
anova.p <- round(anova(mod1)[,5][1], digits=3)
boxplot(map.data[[contrast.name]][map.data$chr=='Chr01'], map.data[[contrast.name]][map.data$chr!='Chr01'], names=list("Chr01", "Autosome"), col=c(2,3), legend=(paste('anova p =', anova.p)), main =paste('Log 2 ratio', group1, 'to', group2, 'expression'), ylab='Log2 expression ratio', cex.main=1.4, cex.lab=2, cex.axis=1.8,ylim=c(-1.5,1.5))

mod1 <- lm(map.data$ratio ~ map.data$chr=='Chr01')
summary(mod1)
anova.p <- round(anova(mod1)[,5][1], digits=3)

cor.sex <- round(cor(group1.expr[map.data$chr=='Chr01'], group2.expr[map.data$chr=='Chr01']), digits=2)
cor.aut <- round(cor(group1.expr[map.data$chr!='Chr01'], group2.expr[map.data$chr!='Chr01']), digits=2)
plot(group1.expr, group2.expr, type="n", xlab=paste(group1,'expression'), ylab=paste(group2,'expression'), main =('Correlation plot'), cex.main=1.4, cex.lab=1.2,xlim=c(0,15),ylim=c(0,15))
points(group1.expr[map.data$chr=='Chr01'], group2.expr[map.data$chr=='Chr01'], pch=20, col=col2)
points(group1.expr[map.data$chr!='Chr01'], group2.expr[map.data$chr!='Chr01'], pch=20, col=col1)
legend('topleft', inset=0.05, legend=c(paste('Chr 1 = ', cor.sex), paste('Autosomes = ',cor.aut)), pch = 20, col=c(col2, col1), cex=0.8 )
legend('bottomright', inset=0.05, legend=(paste('anova p =', anova.p)), cex=0.8)

#codes for calcualte ratio of chr01 and autosome
chr01 <-map.data$ratio[map.data$chr=='Chr01']
chr02 <-map.data$ratio[map.data$chr=='Chr02']
sexchr <-c(chr01,chr02)
chr03 <-map.data$ratio[map.data$chr=='Chr03']
chr04 <-map.data$ratio[map.data$chr=='Chr04']
chr05 <-map.data$ratio[map.data$chr=='Chr05']
chr06 <-map.data$ratio[map.data$chr=='Chr06']
chr07 <-map.data$ratio[map.data$chr=='Chr07']
chr08 <-map.data$ratio[map.data$chr=='Chr08']
chr09 <-map.data$ratio[map.data$chr=='Chr09']
chr10 <-map.data$ratio[map.data$chr=='Chr10']
auto <-c(chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10)

ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(0,2)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) #W = 1160300, p-value = 0.4
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],auto,exact = FALSE) #W = 846290, p-value = 0.6466
wilcox.test(sexchr,auto,exact = FALSE) #W = 1818200, p-value = 0.828

ggsave(file="expr_ratio_males.pdf", path = "/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/")

#boxplot
boxplot(map.data[[contrast.name]][map.data$chr=='Chr01'], map.data[[contrast.name]][map.data$chr!='Chr01'], names=list("Chr01", "Autosome"), col=c(2,3), legend=(paste('anova p =', anova.p)), main =paste('Log 2 ratio', group1, 'to', group2, 'expression'), ylab='Log2 expression ratio', cex.main=1.4, cex.lab=2, cex.axis=1.8,ylim=c(-1.5,1.5))
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(0,2)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 

#boxplot
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(0,2)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) 

mapLength <- sum(max(map.data$start[map.data$chr=='Chr01']), max(map.data$start[map.data$chr=='Chr02']), max(map.data$start[map.data$chr=='Chr03']), max(map.data$start[map.data$chr=='Chr04']), max(map.data$start[map.data$chr=='Chr05']), max(map.data$start[map.data$chr=='Chr09']), max(map.data$start[map.data$chr=='Chr10']))

chr_1_end <- max(map.data$start[map.data$chr=='Chr01'])
chr_2_end <- chr_1_end + max(map.data$start[map.data$chr=='Chr02'])
chr_3_end <- chr_2_end + max(map.data$start[map.data$chr=='Chr03'])
chr_4_end <- chr_3_end + max(map.data$start[map.data$chr=='Chr04'])
chr_5_end <- chr_4_end + max(map.data$start[map.data$chr=='Chr05'])
chr_6_end <- chr_5_end + max(map.data$start[map.data$chr=='Chr06'])
chr_7_end <- chr_6_end + max(map.data$start[map.data$chr=='Chr07'])
chr_8_end <- chr_7_end + max(map.data$start[map.data$chr=='Chr08'])
chr_9_end <- chr_8_end + max(map.data$start[map.data$chr=='Chr09'])
chr_10_end <- chr_9_end + max(map.data$start[map.data$chr=='Chr10'])

#boxplot for dN, dS
ggplot(map.data, aes(x=chr, y=dN, fill=chr)) +
  scale_fill_manual(values = c("red","grey","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(0,2)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))


plot(c(0, chr_10_end), c(min(map.data$ratio), max(map.data$ratio)), axes=F, lwd=2, xlab="Chromosome", ylab="ratio", type="null", col="black", main=paste(group1, '/', group2, 'ratio'))
axis(2, c(0, chr_10_end, 500,2000,4000))
plot(c(0, chr_10_end), c(0, 2), axes=F, lwd=2, xlab="Chromosome", ylab="ratio", type="null", col="black", main=paste(group1, '/', group2, 'ratio'))
axis(2, c(0, 1, 2,3,4,5))

points(map.data$start[map.data$chr=='Chr01'], map.data$ratio[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_1_end + map.data$start[map.data$chr=='Chr02'], map.data$ratio[map.data$chr=='Chr02'], pch=20, lwd=2, type="p",col="8", main="")
points(chr_2_end + map.data$start[map.data$chr=='Chr03'], map.data$ratio[map.data$chr=='Chr03'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_3_end + map.data$start[map.data$chr=='Chr04'], map.data$ratio[map.data$chr=='Chr04'], pch=20, lwd=2, type="p",col="8", main="")
points(chr_4_end + map.data$start[map.data$chr=='Chr05'], map.data$ratio[map.data$chr=='Chr05'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_5_end + map.data$start[map.data$chr=='Chr06'], map.data$ratio[map.data$chr=='Chr06'], pch=20, lwd=2, type="p",col="8", main="")
points(chr_6_end + map.data$start[map.data$chr=='Chr07'], map.data$ratio[map.data$chr=='Chr07'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_7_end + map.data$start[map.data$chr=='Chr08'], map.data$ratio[map.data$chr=='Chr08'], pch=20, lwd=2, type="p",col="8", main="")
points(chr_8_end + map.data$start[map.data$chr=='Chr09'], map.data$ratio[map.data$chr=='Chr09'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_9_end + map.data$start[map.data$chr=='Chr10'], map.data$ratio[map.data$chr=='Chr10'], pch=20, lwd=2, type="p",col="8", main="")
abline(h=1, col=2)
dev.off()

plot(c(0, chr_1_end), c(0, 3), axes=F, lwd=2, xlab="Chromosome1", ylab="ratio", type="null", col="black", main=paste(group1, '/', group2, 'ratio'))
axis(2, c(0, 1, 2,3))
points(map.data$start[map.data$chr=='Chr01'], map.data$ratio[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="1", main="")
points(chr_1_end + map.data$start[map.data$chr=='Chr01'], map.data$ratio[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="8", main="")
abline(v=116443467, col=2)
abline(h=1,col=1)
#dmrt2 transcript:TRINITY_DN112725_c0_g1_i1,location 116306551.

#calculate means by sliding windows
install.packages("RcppRoll")
library("RcppRoll")
install.packages("zoo")
library("zoo")
install.packages("microbenchmark")
library("microbenchmark")
require(zoo)
map.data <- read.table("/Users/Wen-Juan/git/rana_transcriptome/input/sexchr_transc_finallist_dnds.txt", header=T)
str(map.data)

#pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/exp_ratio_chr02.pdf")
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/exp_ratio_chr01.pdf")
chr1_1_start <-map.data$start[map.data$chr=='Chr01']
chr_1_end <- max(map.data$start[map.data$chr=='Chr01'])
#chr2_start <-map.data$start[map.data$chr=='Chr02']
#chr2_end <- max(map.data$start[map.data$chr=='Chr02'])
#width <- rep(seq(from = 20, to =80, by = 100))
zoo.dat <- zoo(map.data$ratio[map.data$chr=='Chr01'], c(chr1_1_start,chr_1_end))
#zoo.dat <- zoo(map.data$ratio[map.data$chr=='Chr02'], c(chr2_start,chr2_end))
#zoo.dat <- zoo(map.data$dS[map.data$chr=='Chr01'], c(chr1_1_start,chr_1_end))
y <- rollapply(zoo.dat, 40, FUN = mean, align = 'center', na.rm=TRUE) #by=1, 
# calculate mean every 15 values with non-overlapping by groups of 10
par(mar=c(5,5,4,3)+0.1) 
  
#plot(c(0, chr_1_end), c(0, 2), axes=F, lwd=2, xlab="Start position (bp)", ylab="dS", type = "null", col="gray20", cex.axis=1.5, cex.lab=1.5, cex.main=2)
plot(c(600000, chr_1_end), c(0, 2), axes=F, lwd=2, xlab="Position (bp)", ylab="XY/XX gene expression ratio", type = "null", col="gray50", cex.axis=1.5, cex.lab=1.5, cex.main=2)
#plot(c(6000, chr2_end), c(0, 2), axes=F, lwd=2, xlab="Position (bp)", ylab="XY/XX gene expression ratio", type = "null", col="gray50", cex.axis=1.5, cex.lab=1.5, cex.main=2)
axis(1,c(600000, 30000000,60000000,90000000, 120000000, 150000000, 180000000, 200000000))
axis(2,c(0,0.5,1,1.5,2))
#axis(2, c(0, 0.2, 0.4, 0.6, 0.8, 1))
#points(map.data$start[map.data$chr=='Chr01'], map.data$dN[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="gray50", main="",cex.axis=1.5)
#points(map.data$start[map.data$chr=='Chr01'], map.data$dS[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="gray50", main="")
#points(map.data$start[map.data$chr=='Chr01'], map.data$ratio[map.data$chr=='Chr01'], pch=20, lwd=2, type="p",col="8", main="")
points(map.data$start[map.data$chr=='Chr02'], map.data$ratio[map.data$chr=='Chr02'], pch=20, lwd=2, type="p",col="8", main="")
abline(v=116443467, col="blue",lwd=2,lty=3)
lines(y,col="blue",lwd=4)
dev.off()

