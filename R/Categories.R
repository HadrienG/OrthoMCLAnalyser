#!/usr/bin/env Rscript

library("ggplot2")
library("argparser")


#Args Parser
# Input should be InfosSummary.txt
parser <- arg.parser("MeanCoreVsPan.R")
parser <- add.argument(parser,c("--input", "--output"),
help=c("input directory", "output file"),
short=c("-i", "-o"))
args <- parse.args(parser)
input=args$input
output=args$output


PanCog = read.table(paste(input,"/PanCog.tmp.txt", sep=""), sep="\t")
CoreCog = read.table(paste(input,"/CoreCog.tmp.txt", sep=""), sep="\t")


PanFreq = as.data.frame(table(PanCog$V2))
PanFreq$Freq = prop.table(PanFreq$Freq)
CoreFreq = as.data.frame(table(CoreCog$V2))
CoreFreq$Freq = prop.table(CoreFreq$Freq)
PanFreq$Quali = "PanGenome"
CoreFreq$Quali = "CoreGenome"

CellularProcess = c("D","M","N","O","T","U","V","W","Y","Z")
Storage = c("A","B","J","K","L")
Metabolism = c("C","E","F","G","H","I","P","Q")
Poorly = c("R","S")

Cog.Categories = c("Cellular processes and signaling", "Information storage and processing", "Metabolism", "Poorly characterized")
Freq = c(sum(PanFreq[ PanFreq$Var1 %in% CellularProcess,]$Freq), sum(PanFreq[ PanFreq$Var1 %in% Storage,]$Freq),
sum(PanFreq[ PanFreq$Var1 %in% Metabolism,]$Freq), sum(PanFreq[ PanFreq$Var1 %in% Poorly,]$Freq))

Pan.cat.data = data.frame(Cog.Categories, Freq)

Freq = c(sum(CoreFreq[ CoreFreq$Var1 %in% CellularProcess,]$Freq), sum(CoreFreq[ CoreFreq$Var1 %in% Storage,]$Freq),
sum(CoreFreq[ CoreFreq$Var1 %in% Metabolism,]$Freq), sum(CoreFreq[ CoreFreq$Var1 %in% Poorly,]$Freq))

Core.cat.data = data.frame(Cog.Categories, Freq)
Core.cat.data$quali = "Core"
Pan.cat.data$quali = "Pan"
Pan.cat.data$Freq = prop.table(Pan.cat.data$Freq)
Core.cat.data$Freq = prop.table(Core.cat.data$Freq)

total <- rbind(Core.cat.data, Pan.cat.data)

Graph = paste(output,"/Categories.pdf", sep="")
pdf(Graph)
ggplot(data = total, aes(x = quali, y = Freq, fill = Cog.Categories)) +
  geom_bar(stat="identity") + scale_fill_brewer(palette = "Set2") +
  labs(x="", y="Coding Sequences", fill="") +
  ggtitle(bquote(atop(.("Functional annotations \nof Core and Pan genome","")))) +
  theme(legend.position="bottom", legend.direction="vertical")
invisible(dev.off())

COLORLIST = c("#EFCCE9", "#C8EC87", "#82DDC9", "#EBB38D", "#E2CA68", "#BFCAB6",
"#F3A9A8", "#B3C7EF", "#95C2CB", "#97E4AE", "#D2F0AC", "#E3AD68", "#B9BE72",
"#D4E2EB", "#E1D6AB", "#CFB5A1", "#E9ABCA", "#F59A8A", "#E7D48F", "#9EC99D",
"#9CE5EA", "#A6DB8C", "#C4B9CE", "#DEB2B4", "#D9EED3", "#E0EA84")

Graph = paste(output,"/Subcategories.pdf", sep="")
pdf(Graph)
ggplot(data = total, aes(x = Quali, y = Freq, fill = Var1)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = COLORLIST) +
  labs(x="", y="Coding Sequences", fill="")
invisible(dev.off())
