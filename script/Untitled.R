

violinData <- data.frame(subset(restab_logCPM_subset, grepl("Chr", chr))["chr"], subset(restab_logCPM_subset, grepl("Chr", chr))[paste('logFC.',colnames(cmat)[k], sep="")])
p = ggplot(violinData, aes_string(factor(violinData$chr),paste('logFC.',colnames(cmat)[k], sep="")), ylim)
p + geom_violin(scale="width", fill=16, trim=F) + geom_boxplot(width=0.15, col=4) + geom_hline(yintercept=-1,lty=3) + geom_hline(yintercept=1, lty=3) + geom_hline(yintercept=2, lty=2)+ geom_hline(yintercept=-2, lty=2) +  scale_y_continuous(breaks=c(-10, -5, -2, -1, 0, 1, 2, 5, 10)) + labs(title=paste('logFC.',colnames(cmat)[k],' DE', sep=""), x='Chromosome', y='logFC') + theme(panel.background = element_rect(fill='white'))  + theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"), axis.title.y = element_text(colour = "black"))
ggsave(file.path(outpath, paste('violin_per_chr_DE_', colnames(cmat)[k], '_', sub_analyse, '.pdf', sep="")), width=8, height=8)
