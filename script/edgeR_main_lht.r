# Load packages, setup colours
library(edgeR)
library(gplots)
library(ggplot2)
library(dynamicTreeCut)
library(statmod)

cbred <- 2 # '#D55E00'
cbblue <- '#0072B2'
cborange <- '#E69F00'
cbgreen <- '#009E73'
cbpink <- '#CC79A7'
cblblue <- '#56B4E9'
cdyellow <- '#F0E442'

col1 <- rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
col1b <- rgb(red = 0, green = 0, blue = 0, alpha = 0.3)
col2 <- rgb(red = 1, green = 0, blue = 0, alpha = 0.6)
col3 <- rgb(red = 0, green = 1, blue = 0, alpha = 0.6)
col4 <- rgb(red = 0, green = 0, blue = 1, alpha = 0.6)
col5 <- rgb(red = 0.7, green = 0, blue = 0.5, alpha = 0.6)

# Setup input (default is command line)
args <- commandArgs(trailingOnly = TRUE)
sub_analyse = paste(args[1])
FDR2use = as.numeric(paste(args[2]))

# example
# sub_analyse <- 'Amm27'
# FDR2use  <- 0.05

datapath <- "~/my_postdoc/useful_scripts/Rana_Transcriptome/input/"
outpath <- paste("~/my_postdoc/useful_scripts/Rana_Transcriptome/output/", sub_analyse, sep="")
dir.create(file.path(outpath))

annotation <- read.delim(file.path(datapath, "Am_annotation.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE) # BRM_annotation.txt annotation_out.txt

count <- read.table(file.path(datapath, paste(sub_analyse,'_count.txt', sep="")), header=T, row.names=1)
count <- round(count, digits=0)
design <- read.table(file.path(datapath, paste(sub_analyse,'_design.txt', sep="")), header=T)

model.formula <- as.formula("~0+group")
dmat <- model.matrix(model.formula,data=as.data.frame(design))
dgl <- DGEList(counts=count, group=design$group, genes=annotation)
paste("all transcripts:", nrow(dgl))

##filter effects
filter_file <- file.path(paste(outpath, '/',sub_analyse, 'filtering_info.txt', sep=""))
ave0 <- dgl[aveLogCPM(dgl) >= 0,]
ave1 <- dgl[aveLogCPM(dgl) >= 1,]
sum2 <- dgl[rowSums(cpm(dgl)) >= 2,]
sum10 <- dgl[rowSums(cpm(dgl)) >= 10,]

write(paste("Comp\tRemaining\tMin\t1Q\tMedian\tMean\t3Q\tMax"), filter_file)
write(paste("All_"), filter_file, append=T)
write(paste(nrow(dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(dgl)/ncol(dgl))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave0_"), filter_file, append=T)
write(paste(nrow(ave0), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave0)/ncol(ave0))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Ave1_"), filter_file, append=T)
write(paste(nrow(ave1), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(ave1)/ncol(ave1))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum2_"), filter_file, append=T)
write(paste(nrow(sum2), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum2)/ncol(sum2))), filter_file, append=T, sep='\t', ncol=6)
write(paste("Sum10_"), filter_file, append=T)
write(paste(nrow(sum10), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(sum10)/ncol(sum10))), filter_file, append=T, sep='\t', ncol=6)

# compile nice table with perl -pe 's/_\n/\t/g' #_filtering_info.txt

# spec <- dgl[rowSums(cpm(dgl)>=2) > 3,] # The gene must be expressed in at least 3 libaries (remove sex specific genes)
dgl <- dgl[aveLogCPM(dgl) > 0,] # filter by average reads
#dgl <- dgl[rowSums(cpm(dgl)>=2) > 3,] #the 5 means there are 10 libraries in total with 5 males and 5 females
dgl <- dgl[rowSums(cpm(dgl)>1) >= 3,]

write(paste("dgl"), filter_file, append=T)
write(paste(nrow(dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(dgl)/ncol(dgl))), filter_file, append=T, sep='\t', ncol=6)

summary(aveLogCPM(dgl))
summary(rowSums(dgl$count))

# MDS ggplot
## colours
col.M23   <- rgb(206/255, 34/255, 43/255, 3/4)
col.F23  <- rgb(206/255, 34/255, 43/255, 3/4)
#col.U23  <- rgb(206/255, 34/255, 43/255, 3/4)
col.M27 <- rgb(209/255, 127/255, 21/255, 3/4)
col.F27 <- rgb(209/255, 127/255, 21/255, 3/4)
#col.U27 <- rgb(209/255, 127/255, 21/255, 3/4)
col.M31 <- rgb(23/255, 87/255, 120/255, 3/4)
col.F31 <- rgb(23/255, 87/255, 120/255, 3/4)
#col.U31 <- rgb(23/255, 87/255, 120/255, 3/4)
col.M43 <- rgb(88/255, 135/255, 37/255, 3/4)
col.F43 <- rgb(88/255, 135/255, 37/255, 3/4)
#col.R43 <- rgb(88/255, 135/255, 37/255, 3/4) #for tvedora
col.M46 <- rgb(113/255, 250/255, 241/255, 3/4)
col.F46 <- rgb(113/255, 250/255, 241/255, 3/4)
#col.R46 <- rgb(113/255, 250/255, 241/255, 3/4)

pdf(file.path(outpath,paste('MDS_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
# design$group
y <- dgl
colnames(y) <- paste(colnames(y), design$group, sep="\n")
#cols = c(col.M23, col.F23,col.M27,col.F27, col.M31, col.F31,col.M43,col.F43,col.R43,col.M46,col.F46,col.R46)#for tvedora
#cols = c(col.M23, col.F23,col.U23,col.M27,col.F27,col.U27, col.M31, col.F31,col.U31,col.M43,col.F43,col.M46,col.F46,col.MF46) #for argovie
cols = c(col.M23, col.F23,col.M27,col.F27, col.M31, col.F31,col.M43,col.F43,col.M46,col.F46) #for ammarnas
#pchs = c(18,5,18,5,18,5,18,5,16,18,5,16) #for tvedora
pchs = c(18,5,18,5,18,5,18,5,18,5)
plotMDS(y, pch=pchs[design$group],col=cols[design$group], cex=2.0, main="Ammarnas MDS plot",cex.main=1, cex.lab=1,lty=2, lwd=3)
legend('bottom', inset=0.02, legend=levels(design$group), pch = pchs, col=cols )
plotMDS(y, cex=1, col=cols[design$group], main=paste(sub_analyse,"MDS plot"))
plotMDS(y, cex=1, pch=as.numeric(y$samples$group), col=as.numeric(y$samples$group), main=paste(sub_analyse,"MDS plot"))
legend('bottom', inset=0.02, legend=levels(design$group), pch = as.numeric(y$samples$group), col=as.numeric(y$samples$group) )
dev.off()

# estimate data normalisation factors and dispersion
xcpm <- mglmOneGroup(dgl$counts)     # computing a logCPM for making dispersion plot
dgl <- calcNormFactors(dgl)
dgl <- estimateGLMCommonDisp(dgl, dmat)
dgl <- estimateGLMTrendedDisp(dgl, dmat, min.n=1000)
dgl <- estimateGLMTagwiseDisp(dgl, dmat)

## dispersion plot
pdf(file.path(outpath, paste('Dispersion_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste(sub_analyse," dispersion", sep=""))
if(!is.null(dgl$trended.dispersion)) points(xcpm,dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=dgl$common.dispersion ,col=cblblue, lwd=2)
legend("topright", c("Common","Trended","Tagwise"), pch=16, col=c(cblblue, cbgreen, "black"), title="Dispersion")
dev.off()

##  fit the data model ------------------------------------------
#fitres <- glmFit(dgl, dmat) #old model fit
#fitres <- glmQLFit(dgl, dmat, robust = TRUE)
fitres <- glmFit(dgl, dmat, robust = TRUE)

x <- read.delim(paste(datapath, sub_analyse, "_matrix.txt", sep=""), sep="\t", header=T)
sortedX <- data.frame(x[order(x$model_coefficients, decreasing=F),])
cmat <- as.matrix(sortedX[,-1])
colnames(cmat)[1] <- colnames(sortedX[2])
rownames(cmat) <- as.character(sortedX[,1])
# cmat

#Contrast fit and test results
lrtres <- list()
#for(k in 1:ncol(cmat)) lrtres[[k]] <- glmQLFTest(fitres, contrast=cmat[,k])
for(k in 1:ncol(cmat)) lrtres[[k]] <- glmLRT(fitres, contrast=cmat[,k])

logFC <- NULL
PV <- NULL
FDR <- NULL
for(k in 1:ncol(cmat)) {PV <- cbind(PV, lrtres[[k]]$table[,"PValue"])
FDR= cbind(FDR, p.adjust(PV[,k],method="BH"))
logFC= cbind(logFC, lrtres[[k]]$table[,"logFC"])
}

xcpm <- lrtres[[1]]$table[,"logCPM"]
allzeros <- which(rowSums(cpm(dgl)) < 3,)
allused <- which(rowSums(cpm(dgl)) >= 3,)

cname <-colnames(cmat)
colnames(logFC) <- paste("logFC", cname,sep=".")
colnames(PV) <- paste("PV", cname, sep=".")
colnames(FDR) <- paste("FDR", cname, sep=".")
rownames(logFC) <- rownames(PV) <- rownames(FDR) <- rownames(fitres$coefficients)

## Make results table
idxzeros <- allzeros

restab <- data.frame(rownames(dgl$counts),dgl$genes[,2], dgl$genes[,3], dgl$genes[,4], dgl$genes[,5],logCPM=xcpm,logFC,PV,FDR) # note start and end are the same in annotation file
#if no go term, no the fifth column
#restab <- data.frame(rownames(dgl$counts),dgl$genes[,2], dgl$genes[,3], dgl$genes[,4], logCPM=xcpm,logFC,PV,FDR) # note start and end are the same in annotation file
colnames(restab)[1] <- 'gid'
colnames(restab)[2] <- 'gname'
colnames(restab)[3] <- 'chr'
colnames(restab)[4] <- 'start'
colnames(restab)[5] <- 'end'

write.table(restab$gid, file=file.path(outpath, paste(sub_analyse, "_used_gene_names.txt", sep="")), quote=F, row.names=F, sep="\t")

## pvalue histogram
pdf(file.path(outpath, paste('p_hist_',FDR2use, '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
npanel <- ncol(logFC)
np <- ceiling(sqrt(npanel))
if(np*(np-1)>= npanel) mfcol <- c(np-1,np) else
{if((np-1)*(np+1)>= npanel) mfcol <- c(np+1,np-1) else mfcol <- c(np,np)}
par(mfcol=mfcol)
for(k in 1:npanel) {
  hist(PV[,k],n=100,xlab="P-value",main=colnames(logFC)[k])
}
dev.off()

## whinin group pairwise scatter plot
wx <- dgl$counts
wg <- as.character(dgl$samples$group)
wug <- unique(wg)
wn <- length(wug)
for (k in 1:wn) {
  ix <- wg %in% wug[k]
  xmat <- log2(wx[,ix])
  pdf(file.path(outpath, paste('pairwise_raw_count_', wug[k], '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
  if (sum(ix) > 1 ) {
    pairs(xmat,pch=16,cex=0.4,main=wug[k])
  }
  #	dev.copy(pdf,file.path(outpath,'courtship_out',paste(sub_analyse,wug[k],'_pairwise_raw_count.pdf', sep="")), width=8, height=8)
  dev.off()
}


##Export number of DE genes table
restab_frame <- as.data.frame(restab)

for(logFC_use in c(3, 2, 1, 0) ) {
  de.yes.no <- FDR < FDR2use & abs(logFC) > logFC_use
  if (ncol(cmat) == 1) {
    de4 <- which((FDR[,1] < FDR2use) == 0 & abs(logFC[,1] > logFC_use) != 0) } else {
      de4 <- which(rowSums(FDR[,c(1:ncol(FDR))] < FDR2use) == 0 & rowSums(abs(logFC[,c(1:ncol(FDR))]) > logFC_use) != 0 )
    }

de.yes.no[de4,] <- FALSE
deidx <- ii <- rowSums(de.yes.no) > 0
delabel <- (sign(logFC)*de.yes.no)[ii,]
delabel[is.na(delabel)] <- 0
combinedp <- 1; for(k in 1:ncol(PV)) combinedp <- combinedp*PV[,k]
deidx <- ii

wtable_1 <- restab[deidx,] # this is the result table of DE genes.
demat <- as.matrix(logFC[deidx,])
write.table(wtable_1, file=file.path(outpath, paste('de_',FDR2use, '_',logFC_use, '_', sub_analyse,'.txt', sep="")), quote=F, row.names=F, sep='\t')

if (nrow(wtable_1) > 1) {
  if (ncol(cmat) == 1) {
    wtable_3 <- rbind(NUM_DE=sum(abs(delabel)), NUM_UP_DE =sum(delabel>0), NUM_DOWN_DE =sum(delabel<0)) } else {
      wtable_3 <- rbind(NUM_DE=colSums(abs(delabel)), NUM_UP_DE =colSums(delabel>0), NUM_DOWN_DE =colSums(delabel<0))
    }
wtable_3 <- cbind(DE_numbers=rownames(wtable_3), wtable_3)
write.table(wtable_3, file=file.path(outpath, paste('Number_de_',FDR2use, '_',logFC_use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
  }
}

#logFC plot
fc <- logFC
fc[(fc)>10] <- 10
fc[ fc < -10] <- -10
par(mfcol=mfcol)

for(k in 1:ncol(cmat)) {
pdf(file.path(outpath, paste('FC-CPMplot_',FDR2use, '_', sub_analyse,'_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=8, height=8)
par(mar=c(5,6,3,2)+0.1)
ylab <- colnames(logFC)[k]
deix <- which(de.yes.no[,k])
maPlot(x=NULL,y=NULL,ylim=c(-10,10),logAbundance= xcpm, logFC = fc[,k], xlab = bquote(paste(log^2, CPM)), ylab = paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' - ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), de.tags= deix, pch = 19, cex = 0.3, smearWidth = 0.5, cex.axis=1.8, cex.lab=2, panel.first = grid(), smooth.scatter = FALSE, lowess = FALSE, na.rm =TRUE, main = paste('LogFC plot ', strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""))
dev.off()
}

# generation of subsets
# colnames(restab_frame)
# tail(restab_frame)
# for(k in 6:ncol(restab_frame)) restab_frame[,k] <- as.numeric(levels(restab_frame[,k]))[restab_frame[,k]]
# str(restab_frame)

## correlation plots, not v useful
# par(mfrow=c(2,2))
# plot(restab[,13], restab[,7], xlab="m vs f", ylab="XXm vs XXf")
# plot(restab[,13], restab[,8], xlab="m vs f", ylab="XYm vs XXf")
# plot(restab[,13], restab[,9], xlab="m vs f", ylab="XYm vs XXm")
# plot(restab[,13], restab[,10], xlab="XXm vs XXf", ylab="XYm vs XXf")
# dev.copy(pdf,file.path(outpath, 'sex bias plot body.pdf'), width=8, height=8)
# dev.off()

# Physical map basics
if (!is.numeric(restab$start)) {restab$start <- as.numeric(levels(restab$start))[restab$start] }
annotation$chr <- as.factor(annotation$chr)
# annotation$start <- annotation$position
annotation$start <- as.numeric(annotation$start)

mapData <- subset(restab, restab$chr!='NA')
mapData2 <- subset(mapData, mapData$start!='NA')
nrow(mapData2)
annotation2 <- subset(annotation, annotation$chr!='NA' & annotation$start!='NA')

mapLength <- sum(max(mapData2$start[mapData2$chr=='Chr01']), max(mapData2$start[mapData2$chr=='Chr02']), max(mapData2$start[mapData2$chr=='Chr03']), max(mapData2$start[mapData2$chr=='Chr04']), max(mapData2$start[mapData2$chr=='Chr05']), max(mapData2$start[mapData2$chr=='Chr09']), max(mapData2$start[mapData2$chr=='Chr10']))

chr_1_end <- max(mapData2$start[mapData2$chr=='Chr01'])
chr_2_end <- chr_1_end + max(mapData2$start[mapData2$chr=='Chr02'])
chr_3_end <- chr_2_end + max(mapData2$start[mapData2$chr=='Chr03'])
chr_4_end <- chr_3_end + max(mapData2$start[mapData2$chr=='Chr04'])
chr_5_end <- chr_4_end + max(mapData2$start[mapData2$chr=='Chr05'])
chr_6_end <- chr_5_end + max(mapData2$start[mapData2$chr=='Chr06'])
chr_7_end <- chr_6_end + max(mapData2$start[mapData2$chr=='Chr07'])
chr_8_end <- chr_7_end + max(mapData2$start[mapData2$chr=='Chr08'])
chr_9_end <- chr_8_end + max(mapData2$start[mapData2$chr=='Chr09'])
chr_10_end <- chr_9_end + max(mapData2$start[mapData2$chr=='Chr10'])


list_de <- list()
list_nonde <- list()

#results for different logFC
for(logFC_use in c(3, 2, 1, 0) ) {
for(k in 1:ncol(cmat)) {
	list_de[[k]] <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use & abs(restab[[paste('logFC.',colnames(cmat)[k], sep="")]]) > logFC_use)
	list_nonde[[k]] <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] > FDR2use | abs(restab[[paste('logFC.',colnames(cmat)[k], sep="")]]) < logFC_use)
}

# Chisq comparison Chr1 with others
for(k in 1:npanel) {
de_chr1 <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr01'] !='na', na.rm=T)
de_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('scaf*', list_de[[k]]$chr) & list_de[[k]]$chr!='Chr01'] !='na', na.rm=T)
nonde_chr1 <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr01'] !='na', na.rm=T)
nonde_chrother <- sum(list_nonde[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('scaf*', list_nonde[[k]]$chr) & list_nonde[[k]]$chr!='Chr01'] !='na', na.rm=T)

up_chr1 <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr01'] > 0, na.rm=T)
up_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('scaf*', list_de[[k]]$chr) & list_de[[k]]$chr!='Chr01'] > 0, na.rm=T)
down_chr1 <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr01'] < 0, na.rm=T)
down_chrother <- sum(list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]][!grepl('scaf*', list_de[[k]]$chr) & list_de[[k]]$chr!='Chr01'] < 0, na.rm=T)

total_significant <- de_chr1 + de_chrother
total_chr1 <- sum(de_chr1 + nonde_chr1)
total_chrother <- sum(de_chrother + nonde_chrother)
exp_chr1 <- round(total_significant * total_chr1/(total_chr1+total_chrother))
exp_chrother <- round(total_significant * total_chrother/(total_chr1+total_chrother))

if (sum(de_chr1 + de_chrother) > 10 & de_chr1 != 0) { # if 0 error in chisqtest stops the script

Obs <- c(de_chr1, de_chrother)
Exp <- c(exp_chr1, exp_chrother)
chi_obs_data <- rbind(Obs, Exp)
colnames(chi_obs_data) <- c('Chr1', 'Chr other')
chi_test <- chisq.test(c(de_chr1, de_chrother), p=c(total_chr1/(total_chr1+total_chrother), total_chrother/(total_chr1+total_chrother)))
pasteX <- paste('chisq =', round(chi_test$p.value, digits=2), sep=" ")

pdf(file.path(outpath, paste('Chisq_FDR', FDR2use, '_logFC_', logFC_use, sub_analyse, '_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=16, height=6)
par(mfrow=c(1,3))
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
mp <- barplot(chi_obs_data, col=c(cbgreen,cbblue), main = paste('DE ',strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], '\nFDR=',FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.4, cex.lab=1.2, ylab="Gene number")
if (chi_test$p.value < 0.001) {title(xlab=pasteX, font.lab=2)} else {title(xlab=pasteX, font.lab=1)}
vals <- cbind(chi_obs_data[1,1], chi_obs_data[2,1], chi_obs_data[1,2], chi_obs_data[2,2])
names(vals) <- LETTERS[1:4]
text(mp, vals, labels = vals, pos = 1, col="white", cex=1.5)
if (sum(vals[1], vals[2]) > sum(vals[3], vals[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
if (sum(vals[1], vals[2]) < sum(vals[3], vals[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
}

if (sum(up_chr1 + up_chrother) > 10 & up_chr1 != 0) { # if 0 error in chisqtest stops the script
total_significant_up <- up_chr1 + up_chrother
exp_chr1_up <- round(total_significant_up * total_chr1/(total_chr1+total_chrother))
exp_chrother_up <- round(total_significant_up * total_chrother/(total_chr1+total_chrother))
Obs_up <- c(up_chr1, up_chrother)
Exp_up <- c(exp_chr1_up, exp_chrother_up)
chi_obs_data_up <- rbind(Obs_up, Exp_up)
colnames(chi_obs_data_up) <- c('Chr1', 'Chr other')
chi_test_up <- chisq.test(c(up_chr1, up_chrother), p=c(total_chr1/(total_chr1+total_chrother), total_chrother/(total_chr1+total_chrother)))
pasteX_up <- paste('chisq =', round(chi_test_up$p.value, digits=2), sep=" ")
mp_up <- barplot(chi_obs_data_up, col=c(cbgreen,cbblue), main = paste('Up ',strsplit(colnames(cmat)[k], '\\.')[[1]][1], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.8, cex.lab=1, ylab="Gene number")
if (chi_test_up$p.value < 0.001) {title(xlab=pasteX_up, font.lab=2)} else {title(xlab=pasteX_up, font.lab=1)}
vals_up <- cbind(chi_obs_data_up[1,1], chi_obs_data_up[2,1], chi_obs_data_up[1,2], chi_obs_data_up[2,2])
names(vals_up) <- LETTERS[1:4]
if (sum(vals_up[1], vals_up[2]) > sum(vals_up[3], vals_up[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
if (sum(vals_up[1], vals_up[2]) < sum(vals_up[3], vals_up[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
text(mp_up, vals_up, labels = vals_up, pos = 1, col="white", cex=1.5)
}

if (sum(down_chr1 + down_chrother) > 10 & down_chr1 != 0) { # if 0 error in chisqtest stops the script
total_significant_down <- down_chr1 + down_chrother
exp_chr1_down <- round(total_significant_down * total_chr1/(total_chr1+total_chrother))
exp_chrother_down <- round(total_significant_down * total_chrother/(total_chr1+total_chrother))
Obs_down <- c(down_chr1, down_chrother)
Exp_down <- c(exp_chr1_down, exp_chrother_down)
chi_obs_data_down <- rbind(Obs_down, Exp_down)
colnames(chi_obs_data_down) <- c('Chr1', 'Chr other')
chi_test_down <- chisq.test(c(down_chr1, down_chrother), p=c(total_chr1/(total_chr1+total_chrother), total_chrother/(total_chr1+total_chrother)))
pasteX_down <- paste('chisq =', round(chi_test_down$p.value, digits=2), sep=" ")
mp_down <- barplot(chi_obs_data_down, col=c(cbgreen,cbblue), main = paste('Down ',strsplit(colnames(cmat)[k], '\\.')[[1]][2], '\nFDR=', FDR2use, ' logFC=', logFC_use, sep=""), beside=T, cex.main=1.8, cex.lab=1, ylab="Gene number")
if (chi_test_down$p.value < 0.001) {title(xlab=pasteX_down, font.lab=2)} else {title(xlab=pasteX_down, font.lab=1)}
vals_down <- cbind(chi_obs_data_down[1,1], chi_obs_data_down[2,1], chi_obs_data_down[1,2], chi_obs_data_down[2,2])
names(vals_down) <- LETTERS[1:4]
text(mp_down, vals_down, labels = vals_down, pos = 1, col="white", cex=1.5)
if (sum(vals_down[1], vals_down[2]) > sum(vals_down[3], vals_down[4])) {legend('topright', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
if (sum(vals_down[1], vals_down[2]) < sum(vals_down[3], vals_down[4])) {legend('topleft', inset=0.05, legend=rownames(chi_obs_data), pch =15, col=c(cbgreen,cbblue) )}
mtext("Sex chromosome vs autosomes", outer = TRUE, cex = 1.5)
}

if (sum(de_chr1 + de_chrother) > 0) {dev.off()}

# Export DE gene names
write.table(list_de[[k]]$gname, file=file.path(outpath, paste('gname_FDR_',FDR2use, '_logFC', logFC_use, '_', colnames(cmat)[k],'_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
write.table(list_de[[k]]$gene, file=file.path(outpath, paste('gene_FDR_', FDR2use, '_logFC', logFC_use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')


# physical map
par(mar=c(5,5,4,3))
pdf(file.path(outpath, paste('Chr_location_FDR_', FDR2use, '_logFC_', logFC_use, '_',colnames(cmat)[k], '_',  sub_analyse, '.pdf', sep="")), width=12, height=8)
plot(c(0, chr_10_end), c(min(1-log(restab[[paste('FDR.',colnames(cmat)[k], sep="")]])), max(1-log(restab[[paste('FDR.',colnames(cmat)[k], sep="")]]))), axes=F, lwd=2, xlab="Chromosome", ylab="1-log(p)", type="null", col="black", main=paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2]))
axis(2, c(0, chr_10_end, 5,10,15))

points(list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr01'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr01']), pch=20, lwd=2, type="p",col="1", main="")
points(chr_1_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr02'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr02']), pch=20, lwd=2, type="p",col="8", main="")
points(chr_2_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr03'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr03']), pch=20, lwd=2, type="p",col="1", main="")
points(chr_3_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr04'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr04']), pch=20, lwd=2, type="p",col="8", main="")
points(chr_4_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr05'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr05']), pch=20, lwd=2, type="p",col="1", main="")
points(chr_5_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr06'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr06']), pch=20, lwd=2, type="p",col="8", main="")
points(chr_6_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr07'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr07']), pch=20, lwd=2, type="p",col="1", main="")
points(chr_7_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr08'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr08']), pch=20, lwd=2, type="p",col="8", main="")
points(chr_8_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr09'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr09']), pch=20, lwd=2, type="p",col="1", main="")
points(chr_9_end + list_nonde[[k]]$start[list_nonde[[k]]$chr=='Chr10'], 1-log(list_nonde[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_nonde[[k]]$chr=='Chr10']), pch=20, lwd=2, type="p",col="8", main="")

# abline (h = c(300), col='blue', lty='dashed')

chrPosition <- c(chr_1_end/2, chr_1_end/2 + chr_2_end/2, chr_2_end/2 + chr_3_end/2, chr_3_end/2 + chr_4_end/2, chr_4_end/2 + chr_5_end/2, chr_5_end/2 + chr_6_end/2, chr_6_end/2 + chr_7_end/2, chr_7_end/2 + chr_8_end/2, chr_8_end/2 + chr_9_end/2, chr_9_end/2 + chr_10_end/2)

chromosomes<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
axis (side=1, lty=0, at = chrPosition, cex.axis=.5, las=1, labels=chromosomes)

# outliers
points(list_de[[k]]$start[list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, cex=0.8, type="p",col=col4, main="")
points(list_de[[k]]$start[list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1,cex=0.8,  type="p",col=col2, main="")
points(chr_1_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p", cex=0.8, col=col4, main="")
points(chr_1_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_2_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_2_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_3_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_3_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_4_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_4_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_5_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_5_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_6_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_6_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_7_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_7_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_8_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_8_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
points(chr_9_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] > logFC_use]), pch=24, lwd=1, type="p",col=col4, cex=0.8, main="")
points(chr_9_end + list_de[[k]]$start[list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use], 1-log(list_de[[k]][[paste('FDR.',colnames(cmat)[k], sep="")]][list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.', colnames(cmat)[k], sep="")]] < -logFC_use]), pch=25, lwd=1, type="p",col=col2, cex=0.8, main="")
legend('topleft', inset=0.05, legend=c('Up', 'Down'), pch =c(24,25), col=c(col4,col2), cex=0.6 )
dev.off()
}
}

# Percentage by chromosome
## subset by chromosome
se <- function(x,y) sqrt(x * (1-x)/y)
for(k in 1:ncol(cmat)) {
M_F_A_nonDE_Chr1 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr01'))
M_F_A_nonDE_Chr2 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr02'))
M_F_A_nonDE_Chr3 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr03'))
M_F_A_nonDE_Chr4 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr04'))
M_F_A_nonDE_Chr5 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr05'))
M_F_A_nonDE_Chr6 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr06'))
M_F_A_nonDE_Chr7 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr07'))
M_F_A_nonDE_Chr8 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr08'))
M_F_A_nonDE_Chr9 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr09'))
M_F_A_nonDE_Chr10 <- nrow(subset(list_nonde[[k]], list_nonde[[k]]$chr == 'Chr10'))

M_F_A_Chr1All <- nrow(subset(restab, restab$chr=='Chr01'))
M_F_A_Chr2All <- nrow(subset(restab, restab$chr=='Chr02'))
M_F_A_Chr3All <- nrow(subset(restab, restab$chr=='Chr03'))
M_F_A_Chr4All <- nrow(subset(restab, restab$chr=='Chr04'))
M_F_A_Chr5All <- nrow(subset(restab, restab$chr=='Chr05'))
M_F_A_Chr6All <- nrow(subset(restab, restab$chr=='Chr06'))
M_F_A_Chr7All <- nrow(subset(restab, restab$chr=='Chr07'))
M_F_A_Chr8All <- nrow(subset(restab, restab$chr=='Chr08'))
M_F_A_Chr9All <- nrow(subset(restab, restab$chr=='Chr09'))
M_F_A_Chr10All <- nrow(subset(restab, restab$chr=='Chr10'))

chr_all <-cbind( M_F_A_Chr1All, M_F_A_Chr2All, M_F_A_Chr3All, M_F_A_Chr4All, M_F_A_Chr5All, M_F_A_Chr6All, M_F_A_Chr7All, M_F_A_Chr8All, M_F_A_Chr9All, M_F_A_Chr10All)

M_F_A_Chr1up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr1down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr01' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr2up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr2down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr02' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr3up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr3down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr03' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr4up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr4down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr04' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr5up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr5down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr05' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr6up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr6down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr06' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr7up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr7down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr07' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr8up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr8down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr08' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr9up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr9down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr09' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))
M_F_A_Chr10up <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] > 0))
M_F_A_Chr10down <- nrow(subset(list_de[[k]], list_de[[k]]$chr=='Chr10' & list_de[[k]][[paste('logFC.',colnames(cmat)[k], sep="")]] < 0))

pc_up <- cbind(M_F_A_Chr1up/M_F_A_Chr1All, M_F_A_Chr2up/M_F_A_Chr2All, M_F_A_Chr3up/M_F_A_Chr3All, M_F_A_Chr4up/M_F_A_Chr4All, M_F_A_Chr5up/M_F_A_Chr5All, M_F_A_Chr6up/M_F_A_Chr6All, M_F_A_Chr7up/M_F_A_Chr7All, M_F_A_Chr8up/M_F_A_Chr8All, M_F_A_Chr9up/M_F_A_Chr9All, M_F_A_Chr10up/M_F_A_Chr10All)

up <- cbind(M_F_A_Chr1up, M_F_A_Chr2up, M_F_A_Chr3up, M_F_A_Chr4up, M_F_A_Chr5up, M_F_A_Chr6up, M_F_A_Chr7up, M_F_A_Chr8up, M_F_A_Chr9up, M_F_A_Chr10up)

pc_down <- cbind(M_F_A_Chr1down/M_F_A_Chr1All, M_F_A_Chr2down/M_F_A_Chr2All, M_F_A_Chr3down/M_F_A_Chr3All, M_F_A_Chr4down/M_F_A_Chr4All, M_F_A_Chr5down/M_F_A_Chr5All, M_F_A_Chr6down/M_F_A_Chr6All, M_F_A_Chr7down/M_F_A_Chr7All, M_F_A_Chr8down/M_F_A_Chr8All, M_F_A_Chr9down/M_F_A_Chr9All, M_F_A_Chr10down/M_F_A_Chr10All)

down <- cbind(M_F_A_Chr1down, M_F_A_Chr2down, M_F_A_Chr3down, M_F_A_Chr4down, M_F_A_Chr5down, M_F_A_Chr6down, M_F_A_Chr7down, M_F_A_Chr8down, M_F_A_Chr9down, M_F_A_Chr10down)

non_DE <- cbind(M_F_A_nonDE_Chr1, M_F_A_nonDE_Chr2, M_F_A_nonDE_Chr3, M_F_A_nonDE_Chr4, M_F_A_nonDE_Chr5, M_F_A_nonDE_Chr6,	M_F_A_nonDE_Chr7, M_F_A_nonDE_Chr8, M_F_A_nonDE_Chr9, M_F_A_nonDE_Chr10)

se_up <- se(pc_up, chr_all)
se_down <- se(pc_down, chr_all)

pdf(file.path(outpath, paste('logFC_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=12, height=6)
par(mfrow=c(1,2))
par(mar=c(5,5,4,3))
plot(c(0,10), c(0,max(pc_up+se_up, pc_down+se_down)), type='n', xlab="Chromosome", ylab="% of genes", main=paste('Up ',strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), cex.main=1.8, cex.lab=1.3)
points(seq(1:10), pc_up)
points(seq(1:10), pc_up+se_up, cex=0.5)
points(seq(1:10), pc_up-se_up, cex=0.5)
plot(c(0,10), c(0,max(pc_up+se_up, pc_down+se_down)), type='n', xlab="Chromosome", ylab="% of genes", main=paste('Down ',strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), cex.main=1.8, cex.lab=1.3)
points(seq(1:10), pc_down)
points(seq(1:10), pc_down+se_down, cex=0.5)
points(seq(1:10), pc_down-se_down, cex=0.5)
dev.off()
}

## Output subsets FDR05 for GO annotation
### 5% FDR
for(k in 1:ncol(cmat)) {
FDR05 <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
FDR05_unbiased <- subset(restab, restab[[paste('FDR.',colnames(cmat)[k], sep="")]] > FDR2use)
UP <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0)
DOWN <- subset(FDR05, FDR05[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0)

write.table(list_de[[k]]$gid, file=file.path(outpath, paste('gene_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

write.table(data.frame(gene=FDR05_unbiased$gid, logFC=FDR05_unbiased[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=FDR05_unbiased[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('unbiased_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=UP$gid, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=UP[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('UP_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')
write.table(data.frame(gene=DOWN$gid, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]], FDR=DOWN[[paste('FDR.',colnames(cmat)[k], sep="")]]), file=file.path(outpath, paste('DOWN_',FDR2use, '_', colnames(cmat)[k], '_', sub_analyse, '.txt', sep="")), quote=F, sep='\t')

gopath <- file.path(outpath, paste('GO_',FDR2use, '_', colnames(cmat)[k], sep=""))
dir.create(gopath)

write.table(data.frame(gene=FDR05$gid), file=file.path(gopath,'FDR_genes.txt'), quote=F, row.names=F, sep='\t')
write.table(data.frame(gene=UP$gid, logFC=UP[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")
write.table(data.frame(gene=DOWN$gid, logFC=DOWN[[paste('logFC.',colnames(cmat)[k], sep="")]]), file=file.path(gopath, paste(strsplit(colnames(cmat)[k], '\\.')[[1]][2], '_biased.txt', sep="")), quote=F, row.names=F, sep="\t")

GO_pvalues <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]], logFC=restab[[paste('logFC.',colnames(cmat)[k], sep="")]])
write.table(GO_pvalues, file=file.path(gopath, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_UP <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
GO_pvalues_UP[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] < 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
write.table(GO_pvalues_UP, file=file.path(gopath, "GO_pvalues_UP.txt"), quote=F, row.names=F, sep="\t")

GO_pvalues_DOWN <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]])
GO_pvalues_DOWN[,2][restab[[paste('logFC.',colnames(cmat)[k], sep="")]] > 0 & restab[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use] <- 1
write.table(GO_pvalues_DOWN, file=file.path(gopath, "GO_pvalues_DOWN.txt"), quote=F, row.names=F, sep="\t")
}

## Export results and moderated log-counts-per-million
nc <- cpm(dgl, prior.count=2, log=T)
nc2 <- data.frame(row.names(nc), nc)

colnames(nc2) <- c('ID', as.character(design$sample))
restab_logCPM = merge(restab, nc2, by.x="gid", by.y="ID", all=F )
write.table(restab_logCPM, file=file.path(outpath, paste('LogCPM_',FDR2use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')

# look at make specific genes
# male_specific <- subset(restab_logCPM, restab_logCPM$FDR.m.f < 0.01)
# summary(rowSums(cbind(male_specific$F1O, male_specific$F2O, male_specific$F4O))/3 < -1)

## Violin plots
for(k in 1:ncol(cmat)) {
par(mar=c(5,5,4,3))
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

for(k in 1:ncol(cmat)) {

minUse <- min(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")])
chr1_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr01'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr01' & restab_logCPM$chr !='Chr01'])$p.value
chr2_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr02'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr02' & restab_logCPM$chr !='Chr02'])$p.value
chr3_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr03'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr03' & restab_logCPM$chr !='Chr03'])$p.value
chr4_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr04'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr04' & restab_logCPM$chr !='Chr04'])$p.value
chr5_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr05'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr05' & restab_logCPM$chr !='Chr05'])$p.value
chr6_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr06'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr06' & restab_logCPM$chr !='Chr06'])$p.value
chr7_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr07'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr07' & restab_logCPM$chr !='Chr07'])$p.value
chr8_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr08'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr08' & restab_logCPM$chr !='Chr08'])$p.value
chr9_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr09'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr09' & restab_logCPM$chr !='Chr09'])$p.value
chr10_p <- wilcox.test(restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr10'], restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] & restab_logCPM$chr !='Chr10' & restab_logCPM$chr !='Chr10'])$p.value

	pdf(file.path(outpath, paste('logFC_per_chr_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
boxplot(
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr01'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr02'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr03'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr04'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr05'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr06'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr07'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr08'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr09'],
	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")]  & restab_logCPM$chr =='Chr10'],
type="n", xlab="Chromosome", names=seq(1:10), ylab="logFC", main=paste('logFC per chromosome ',strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), cex.main=1.5, cex.lab=1.3)

if (chr1_p < 0.001) {text(1, minUse, labels='***', cex=1.5)} else if (chr1_p < 0.01) {text(1, minUse, labels='**', cex=1.5)} else if (chr1_p < 0.05) {text(1, minUse, labels='*', cex=1.5)}
if (chr2_p < 0.001) {text(2, minUse, labels='***', cex=1.5)} else if (chr2_p < 0.01) {text(2, minUse, labels='**', cex=1.5)} else if (chr2_p < 0.05) {text(2, minUse, labels='*', cex=1.5)}
if (chr3_p < 0.001) {text(3, minUse, labels='***', cex=1.5)} else if (chr3_p < 0.01) {text(3, minUse, labels='**', cex=1.5)} else if (chr3_p < 0.05) {text(3, minUse, labels='*', cex=1.5)}
if (chr4_p < 0.001) {text(4, minUse, labels='***', cex=1.5)} else if (chr4_p < 0.01) {text(4, minUse, labels='**', cex=1.5)} else if (chr4_p < 0.05) {text(4, minUse, labels='*', cex=1.5)}
if (chr5_p < 0.001) {text(5, minUse, labels='***', cex=1.5)} else if (chr5_p < 0.01) {text(5, minUse, labels='**', cex=1.5)} else if (chr5_p < 0.05) {text(5, minUse, labels='*', cex=1.5)}
if (chr6_p < 0.001) {text(6, minUse, labels='***', cex=1.5)} else if (chr6_p < 0.01) {text(6, minUse, labels='**', cex=1.5)} else if (chr6_p < 0.05) {text(6, minUse, labels='*', cex=1.5)}
if (chr7_p < 0.001) {text(7, minUse, labels='***', cex=1.5)} else if (chr7_p < 0.01) {text(7, minUse, labels='**', cex=1.5)} else if (chr7_p < 0.05) {text(7, minUse, labels='*', cex=1.5)}
if (chr8_p < 0.001) {text(8, minUse, labels='***', cex=1.5)} else if (chr8_p < 0.01) {text(8, minUse, labels='**', cex=1.5)} else if (chr8_p < 0.05) {text(8, minUse, labels='*', cex=1.5)}
if (chr9_p < 0.001) {text(9, minUse, labels='***', cex=1.5)} else if (chr9_p < 0.01) {text(9, minUse, labels='**', cex=1.5)} else if (chr9_p < 0.05) {text(9, minUse, labels='*', cex=1.5)}
if (chr10_p < 0.001) {text(10, minUse, labels='***', cex=1.5)} else if (chr10_p < 0.01) {text(10, minUse, labels='**', cex=1.5)} else if (chr10_p < 0.05) {text(10, minUse, labels='*', cex=1.5)}

	#	restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")][restab_logCPM[paste('logFC.',colnames(cmat)[k], sep="")] > 0 & restab_logCPM$chr !='Chr01' & restab_logCPM$chr !='Chr01'], # autosomes
dev.off()
}

# Heatmaps
colnames(restab_logCPM)
tail(restab_logCPM)
cmat <- data.frame(cmat)

for(k in 1:ncol(cmat)) {
cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]]!=0)
design_subset <- design[design$group %in% row.names(cmat_subset),]

length(as.character(design_subset$sample))
rownames(restab_logCPM) <- restab_logCPM$gid
DE_counts <- subset(restab_logCPM, restab_logCPM[[paste('FDR.',colnames(cmat)[k], sep="")]] < FDR2use)
DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))
# row.names(DE_counts_relevant) <- DE_counts$gname
colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group)
if (nrow(DE_counts_relevant) > 2) {

pdf(file.path(outpath, paste('Heatmap_', FDR2use, '_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
heatmap.2(as.matrix(DE_counts_relevant), col=colorpanel(100, cbred,'white', cbblue), scale="row", key=T, keysize=1.3, density.info="density", trace="none", cexCol=0.5, cexRow=0.6, main=paste(sub_analyse, strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use), srtCol=0, key.title=NA)

cmethods <- c('average', 'ward.D') # complete ward.D average median ward.D2 mcquitty centroid single

for(m in 1:length(cmethods)) {
d <- as.matrix(DE_counts_relevant)
myheatcol <- colorpanel(100, cbred,'white', cbblue) # choose a color palette for the heat map
distmatrix <- as.dist(1-cor(t(d), method="pearson"))
hr <- hclust(distmatrix, method=cmethods[m])  # plot(hr)
mycl <- cutreeDynamic(hr, distM=as.matrix(distmatrix)) # dynamic tree cut
clusterCols <- rainbow(length(unique(mycl))) # get a color palette equal to the number of clusters
 myClusterSideBar <- clusterCols[mycl] # create vector of colors for side bar, old code
myClusterSideBar <- as.character(as.numeric(mycl))
heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.5, cexRow=0.6, RowSideColors = myClusterSideBar, main=paste(sub_analyse, cmethods[m], strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use), srtCol=45, key.title=NA)
legend(x=0.15,y=1.12, legend = unique(mycl), col = unique(as.numeric(mycl)), lty= 1, lwd = 3, cex=.5, title="clusters", xpd=T)
DE_counts[[ cmethods[m] ]] = as.factor(mycl)
gopath_method <- file.path(outpath, paste('GO_',FDR2use, '_', colnames(cmat)[k],'/',cmethods[m], sep=""))
dir.create(gopath_method)

#for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
#FDR05 <- data.frame(gene=DE_counts$gid, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts[[paste('FDR.',colnames(cmat)[k], sep="")]] )
#FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1
#write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')

for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
FDR05 <- data.frame(gene=DE_counts$gid, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts[[paste('FDR.',colnames(cmat)[k], sep="")]] )
FDR05[,3][FDR05[,3] < FDR2use & FDR05[,2] != l] <- 1
GO_pvalues <- data.frame(gene=as.character(restab$gid), FDR=restab[[paste('FDR.',colnames(cmat)[k], sep="")]], logFC=restab[[paste('logFC.',colnames(cmat)[k], sep="")]])
merged1 <- merge(GO_pvalues, FDR05, by.x="gene", by.y="gene", all=T )
merged1$FDR.y[is.na(merged1$FDR.y)] <- as.character(merged1$FDR.x[is.na(merged1$FDR.y)])
write.table(data.frame(gene=merged1[,1], FDR=merged1[,5]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')
  # write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')

}
}
dev.off()
##print(paste(colnames(cmat)[k], 'groups:  median', length(levels(DE_counts$median)), '  ward', length(levels(DE_counts$ward)), '  complete', length(levels(DE_counts$complete)), '  average', length(levels(DE_counts$average)), '  ward2', length(levels(DE_counts$ward.D2)), '  mcquitty', length(levels(DE_counts$mcquitty)), '  centroid', length(levels(DE_counts$centroid)), '  centroid', length(levels(DE_counts$centroid)),sep=" "))

write.table(DE_counts, file=file.path(outpath, paste('clusters',FDR2use, '_', sub_analyse, '_',colnames(cmat)[k], '.txt', sep="")), quote=F, row.names=F, sep='\t')
}
}


# heatmaps all
#for(k in 1:ncol(cmat)) {
#cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]]!=0)
#design_subset <- design[sort(design$group) %in% sort(row.names(cmat_subset)),]

#length(as.character(design_subset$sample))
#rownames(restab_logCPM) <- restab_logCPM$gene
#DE_counts <- restab_logCPM
#DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))
#row.names(DE_counts_relevant) <- NULL

#colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group)
 #if (nrow(DE_counts_relevant) > 2) {

#pdf(file.path(outpath, paste('Heatmap_All_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
#cmethods <- c('median', 'average', 'ward.D', 'complete') # complete ward.D average median ward.D2 mcquitty centroid single

# for(m in 1:length(cmethods)) {
#d <- as.matrix(DE_counts_relevant)
#myheatcol <- colorpanel(100, cbred,'white', cbblue) # choose a color palette for the heat map
#distmatrix <- as.dist(1-cor(t(d), method="pearson"))
#hr <- hclust(distmatrix, method=cmethods[m])  # plot(hr)
#mycl <- cutreeDynamic(hr, distM=as.matrix(distmatrix)) # dynamic tree cut
#clusterCols <- rainbow(length(unique(mycl))) # get a color palette equal to the number of clusters
# # myClusterSideBar <- clusterCols[mycl] # create vector of colors for side bar, old code
#myClusterSideBar <- as.character(as.numeric(mycl))
#heatmap.2(d, col=myheatcol, Rowv=reorder(as.dendrogram(hr), wts=mycl), keysize=1.3, scale="row", density.info="density", trace="none", cexCol=0.7, cexRow=0.6, RowSideColors = myClusterSideBar, main=paste(sub_analyse, cmethods[m], strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], FDR2use), srtCol=45, key.title=NA)
#legend(x=0.15,y=1.12, legend = unique(mycl), col = unique(as.numeric(mycl)), lty= 1, lwd = 3, cex=.5, title="clusters", xpd=T)
#DE_counts[[ cmethods[m] ]] = as.factor(mycl)
#gopath_method <- file.path(outpath, paste('GO_All_', colnames(cmat)[k],'/',cmethods[m], sep=""))
#dir.create(gopath_method)
#
 #for(l in 1:length( levels (DE_counts[[cmethods[m]]] ) ) ) {
#FDR05 <- data.frame(gene=DE_counts$gene, cluster=DE_counts[[cmethods[m]]], FDR=DE_counts[[paste('FDR.',colnames(cmat)[k], sep="")]] )
# write.table(data.frame(gene=FDR05[,1], FDR=FDR05[,3]), file=file.path(gopath_method, paste('cluster_',l ,'_',cmethods[m],'.txt', sep="")), quote=F, row.names=F, sep='\t')
 #}
 #}
 #dev.off()

# print(paste(colnames(cmat)[k], 'groups:  median', length(levels(DE_counts$median)), '  ward', length(levels(DE_counts$ward)), '  complete', length(levels(DE_counts$complete)), '  average', length(levels(DE_counts$average)), '  ward2', length(levels(DE_counts$ward.D2)), '  mcquitty', length(levels(DE_counts$mcquitty)), '  centroid', length(levels(DE_counts$centroid)), '  centroid', length(levels(DE_counts$centroid)),sep=" "))

# write.table(DE_counts, file=file.path(outpath, paste('clusters_All_', sub_analyse, '_',colnames(cmat)[k], '.txt', sep="")), quote=F, row.names=F, sep='\t')
# }
#}


# all - takes long and not too useful so commented out
 #for(k in 1:ncol(cmat)) {
 #pdf(file.path(outpath, paste('Heatmap_all_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
 #cmat_subset <- subset(cmat, cmat[[colnames(cmat)[k]]]!=0)
 #row.names(cmat_subset)
# # str(design)
# design_subset <- design[design$group %in% row.names(cmat_subset),]

# length(as.character(design_subset$sample))

 #DE_counts <- restab_logCPM
 #DE_counts_relevant <- subset(DE_counts, select=as.character(design_subset$sample))
 ## row.names(DE_counts_relevant) <- NULL
# colnames(DE_counts_relevant) <- paste(colnames(DE_counts_relevant), '\n', design_subset$group)
# heatmap.2(as.matrix(DE_counts_relevant), col=colorpanel(100, cbred,'white',cbblue), scale="row", key=T, keysize=1.3, density.info="density", trace="none", cexCol=0.5, cexRow=0.6, main=paste(sub_analyse, strsplit(colnames(cmat)[k], '\\.')[[1]][1], 'vs', strsplit(colnames(cmat)[k], '\\.')[[1]][2], 'All'), srtCol=45, key.title=NA)
 #dev.off()
 #}
