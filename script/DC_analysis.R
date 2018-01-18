## need manual changing
compname <- 'Amg46male' #'tv46' #G46
contrast.name <- 'logFC.XY46.SR46' # logFC.XYm.XXmt
group1 <- 'XY'
group2 <- 'XX'

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)

datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/'
outpath <- ("~/my_postdoc/useful_scripts/Rana_Transcriptome/output/")
kdata <- read.table(file.path(datapath, compname, paste('LogCPM_0.05_', compname, '.txt', sep="")), header=T)
str(kdata)

map.data <- subset(kdata, kdata$start!='NA')
chr1.data <- subset(map.data, map.data$chr=='Chr01') 
chr2.data <- subset(map.data, map.data$chr!='Chr01') 

## need changing, look at design
group1.expr.minus <- (map.data$Am2_463 + map.data$Am4_461 + map.data$Am6_462) / 3 #XY group #G46
group2.expr.minus <- map.data$Am5_461 #XX group #G46


group1.expr <- group1.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- group1.expr/group2.expr

map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

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

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/output/figures/dc_xy_xx_testis.pdf", width=8, height=8)
ggplot(map.data, aes(x=chr, y=ratio, fill=chr)) +
  scale_fill_manual(values = c("red","red","grey","grey","grey","grey","grey","grey","grey","grey")) +
  scale_y_continuous(limits = c(0,2)) + 
  geom_boxplot() +
  labs(x='Chromosome', y='Ratio of XY/XX gene expression') +
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

dev.off()

wilcox.test(map.data$ratio[map.data$chr=='Chr01'],map.data$ratio[map.data$chr!='Chr01'],exact = FALSE) #W = 2004300, p-value = 0.5562
wilcox.test(map.data$ratio[map.data$chr=='Chr02'],auto,exact = FALSE) #W = 1465600, p-value = 0.2975
wilcox.test(sexchr,auto,exact = FALSE) #W = 3147100, p-value = 0.2423


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

