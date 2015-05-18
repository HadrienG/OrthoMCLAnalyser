#!/usr/bin/env Rscript

library("ggplot2")
library("argparser")


#Args Parser
# Input should be OutDir
parser <- arg.parser("MeanCoreVsPan.R")
parser <- add.argument(parser,c("--input", "--output", "--ucore", "--upan"),
help = c("input directory", "output file", "Core not in COG", "Pan not in Cog"),
short = c("-i", "-o", "-c", "-p"))
args <- parse.args(parser)
input = args$input
output = args$output
UCore = strtoi(args$ucore)
UPan = strtoi(args$upan)



PanCog = read.table(paste(input,"/PanCog.tmp.txt", sep=""), sep="\t")
CoreCog = read.table(paste(input,"/CoreCog.tmp.txt", sep=""), sep="\t")


# Create data frames of unknown cogs
V1 <- rep("Unknown", UCore)
V2 <- rep(NA, UCore)
V3 <- rep("Not in COG Database", UCore)
CoreRows = data.frame(V1, V2, V3)

V1 <- rep("Unknown", UPan)
V2 <- rep(NA, UPan)
V3 <- rep("Not in COG Database", UPan)
PanRows = data.frame(V1, V2, V3)

fullCore = rbind(CoreCog, CoreRows)
fullPan = rbind(PanCog, PanRows)

# PanCog --> fullPan
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

COLORLIST = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C",
"#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1",
"#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0",
"#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

PanFreq = as.data.frame(table(fullPan$V3))
PanFreq$Freq = prop.table(PanFreq$Freq)
CoreFreq = as.data.frame(table(fullCore$V3))
CoreFreq$Freq = prop.table(CoreFreq$Freq)
PanFreq$Quali = "PanGenome"
CoreFreq$Quali = "CoreGenome"
total <- rbind(CoreFreq, PanFreq)

Graph = paste(output,"/Subcategories.pdf", sep="")
pdf(Graph)
ggplot(data = total, aes(x = Quali, y = Freq, fill = Var1)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = COLORLIST) +
  labs(x="", y="Coding Sequences", fill="") +
  guides(fill = guide_legend(reverse=TRUE))
invisible(dev.off())
