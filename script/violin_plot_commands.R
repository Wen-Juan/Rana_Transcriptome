data_all <- read.table("~/path/violin_plot_file.txt",stringsAsFactors=F,header=T, row.names=1, sep=",")

# data_all is of the form:
# gene, logFC, group
# where I've divided the data for the two tissues I had into two groups (but you can use it just with one group as well)

p = ggplot(data_all, aes(factor(group),logFC), ylim)
p + geom_violin(scale="width") + geom_hline(yintercept=-1,linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") + theme(panel.background = element_rect(fill='white')) + scale_y_continuous(breaks=c(-10, -5, -1, 0, 1, 5, 10, 15)) +theme(axis.text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(colour = "black"),
         axis.title.y = element_text(colour = "black"))

# I've used geom_hline(yintercept=-1,linetype="dashed") + geom_hline(yintercept=1, linetype="dashed") just to make a mental note of where the male- and female-biased region starts within the violin plot, so I know where to break the color gradient when using Inkscape.

