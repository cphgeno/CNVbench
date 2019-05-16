library(ggplot2)
library(reshape2)

pdf("XXX.pdf", width = 15)
####################
### Prepare data ###
####################
hist <- read.table("hist_data.txt", stringsAsFactors=FALSE, header=TRUE, na.strings = ".")
# Force bin order
hist$bin <- factor(hist$bin, levels = c("100000-Inf","10000-100000","1000-10000","500-1000","0-500"))

plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts, fill=bin)) + facet_grid(~library)
#plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts, fill=library)) + facet_grid(~bin)
#plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts)) + facet_grid(bin~library) + scale_y_log10()
plot + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

for (library in c("WES", "WGS")) {
  lib <- hist[hist$library == library,]
  plot <- ggplot(data = lib) + geom_col(aes(x=sample, y=counts, fill=bin)) + facet_grid(~tool)
  print(plot + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(library))
}





####################
### Prepare data ###
####################
hist <- read.table("summary.txt", stringsAsFactors=FALSE, header=TRUE, na.strings = ".")
# Remove value column
hist <- subset(hist, select=-c(TP,FN,FP,TP_FN,N_truth,precision,recall))
# meld
hist <- melt(hist, value.name="counts")
hist <- plyr::rename(hist, c("variable"="cnv"))

plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts, fill=cnv)) + facet_grid(~library)
#plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts, fill=library)) + facet_grid(~cnv)
#plot <- ggplot(data = hist) + geom_col(aes(x=tool, y=counts)) + facet_grid(cnv~library)
plot + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

for (library in c("WES", "WGS")) {
  lib <- hist[hist$library == library,]
  plot <- ggplot(data = lib) + geom_col(aes(x=sample, y=counts, fill=cnv)) + facet_grid(~tool)
  print(plot + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ggtitle(library))
}

dev.off()

