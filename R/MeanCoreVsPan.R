#!/usr/bin/env Rscript

library("ggplot2")
library("argparser")

#Args Parser
# Input should be InfosSummary.txt
parser <- arg.parser("MeanCoreVsPan.R")
parser <- add.argument(parser,c("--input", "--output"),
help=c("input file", "output file"),
short=c("-i", "-o"))
args <- parse.args(parser)
input=args$input
output=args$output

data=read.table(input,header=T)

meanc = mean(data$CoreSize)
meanp = mean(data$PanSize)

dat = data.frame(
  gen = factor(c("Core","Pan"), levels=c("Core","Pan")),
  size = c(meanc,meanp)
  )

plot.title = ""
plot.subtitle = ""

Graph = paste(output,"/MeanCoreVsPan.pdf", sep="")
pdf(Graph)
ggplot(data=dat, aes(x=gen, y=size, fill=gen)) +
  geom_bar(stat="identity") + ylab("Size (bp)") + xlab("") +
  guides(fill=guide_legend(title="")) +
  ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))))
invisible(dev.off())
