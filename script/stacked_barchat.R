#install and load packages
install.packages("ggplot2")
install.packages("reshape2")

library(ggplot2)
library(reshape2)

#set up working directory
setwd <-'~/Rana_Transcriptome/'
datapath <- '~/Rana_Transcriptome/input/'
results <- '~/Rana_Transcriptome/output/'

#load datasets
sb <- read.table('~/Rana_Transcriptome/input/amm_sexbias.txt', header = TRUE)
str(sb)

sb <- data.frame(with(sb,sb[order(Gosner_stage,Fold_change,Sex_bias),]))
head(sb)

blues <- RColorBrewer::brewer.pal(4, "Blues")
reds <- RColorBrewer::brewer.pal(4, "Reds")

par(mar=c(5,1,3,5))
pdf("~/Rana_Transcriptome/output/figures/sex_bias_alltissues.pdf",height=10, width=12)  

ggplot(data = sb, aes(x=Sex_bias, y= Gene_number, fill=interaction(factor(Fold_change),Sex_bias), label=Gene_number)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(reds,blues), name="Fold change",labels=c("FDR ≤0.05","0< Log2 ≤1","1< Log2 ≤2","2< Log2 ≤3","FDR ≤0.05","0< Log2 ≤1","1< Log2 ≤2","2< Log2 ≤3")) +
  facet_grid(~Gosner_stage) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_discrete(labels=c("F", "M"),name="Sex bias") +
  scale_y_continuous(name = "Number of genes", limits = c(0,6500)) + 
  theme(axis.title.x = element_text(size=15,colour = "black"),axis.title.y = element_text(size=15,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12),axis.text.y = element_text(colour="black",size=12))

dev.off()

