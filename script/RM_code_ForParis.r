##Assigning Transcripts to Chromosomes###
MAChr1 <- subset(mydata, mydata$LG=="1")


## Sort According to Position ####
MAChr1_sort <- MAChr1[order(MAChr1$cf_female),] 


###### Sliding Window Analysis - Judith Visualization ####

library(zoo)

Chr1RM<- rollmean(smooth(MAChr1_sort$WEIGHTED_FST),50)
Chr2RM<- rollmean(smooth(MAChr2_sort$WEIGHTED_FST),50)
Chr3RM<- rollmean(smooth(MAChr3_sort$WEIGHTED_FST),50)
Chr4RM<- rollmean(smooth(MAChr4_sort$WEIGHTED_FST),50)
Chr5RM<- rollmean(smooth(MAChr5_sort$WEIGHTED_FST),50)
Chr6RM<- rollmean(smooth(MAChr6_sort$WEIGHTED_FST),50)
Chr7RM<- rollmean(smooth(MAChr7_sort$WEIGHTED_FST),50)
Chr8RM<- rollmean(smooth(MAChr8_sort$WEIGHTED_FST),50)


testRM <- c(Chr2RM, Chr3RM, Chr4RM, Chr5RM, Chr6RM, Chr7RM, Chr8RM)


myfunction <- function(i){
Info <- sample(i,1,replace=FALSE)
return(Info)
}

my.perm <- c()
for(i in 1:10^3){ my.perm[i] <- myfunction(testRM) }
sorted.perm <- sort(my.perm)
lowCI <- sorted.perm[25]
highCI <- sorted.perm[975]



RMpalette <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")

pdf("/Users/melissatoups/Desktop/Mercurialis/Annua_Prelim_Sex/Exome_Data/For_Genome_Paper/Double.Check.Alignments/by_sliding_window/MA_LG1_WEIGHTED_FST_RM.pdf", width=7,height=5)
MAChr1RM <- rollmean(smooth(MAChr1_sort$cf_female), 50)
MAWEIGHTED_FSTChr1RM <- rollmean(smooth(MAChr1_sort$WEIGHTED_FST),50)
plot(MAChr1_sort$cf_female, MAChr1_sort$WEIGHTED_FST,col=alpha(RMpalette[3], 0.5),pch=20, ylim=c(-0.05,0.3), xlab="Position(CM)", ylab="WEIGHTED_FST",main="Chr1")
lines(MAChr1RM, MAWEIGHTED_FSTChr1RM,type="l",lwd=5, col=RMpalette[5])
abline(h=lowCI,lty=2)
abline(h=highCI,lty=2)
dev.off()

