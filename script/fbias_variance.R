datapath <- '/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/'
outpath <- ("~/my_postdoc/useful_scripts/Rana_Transcriptome/output/Amm_brainissue/")
amm_fbias_anno <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Rana_Transcriptome/input/fbias_xtropannot_all_score.txt",header = T)
str(amm_fbias_anno)
boxplot(score~Group,data=amm_fbias_anno,	xlab="Groups of f-bias", ylab="Blast Score")

group1.expr <- group1.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))
group2.expr <- group2.expr.minus + abs(min(c(group1.expr.minus, group2.expr.minus)))

map.data$ratio <- group1.expr/group2.expr

map.data$ratio[mapply(is.infinite, map.data$ratio)] <- NA

write.table(map.data, "~/my_postdoc/useful_scripts/Rana_Transcriptome/output/map.data", sep="\t")