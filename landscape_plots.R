#!/usr/bin/env Rscript
# Usage: Rscript --vanilla landscape_divergence_plot.R RMSK_landscape_file GenomeSize_summary.txt list_sats figure_name

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("4 arguments must be supplied", call.=FALSE)
} else if (length(args) == 4) {
  filename = args[1]
  sizefile = args[2]
  plotname = args[4]
  list_sats = args[3]
}

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)


aves_bps = data.frame()

# read input files
read.landscape = function(file){
  
  require(data.table)
  data = fread(input = file)
  data = data[,c(1:17)]
  
  return(data)
}
aves = read.landscape(file = filename)

# read genome size summary file (libraries)
genomeSize = fread(sizefile)

# read list of satellite sequences to plot
sats = fread(list_sats, header = FALSE)$V1

# clean dataset
clean.landscape = function(data){
  
  # keep only the columns with: name, size, divergence, ID number
  data = data[,c(10,17,16,15)]
  names(data) = c("Repeat", "Size", "Divergence", "ID")
  data = data[grepl(pattern = "Satellite", x = data$Repeat),]
  elements = strsplit(x = data$Repeat, split = "#")
  data$Subfamily = unlist(sapply(elements, "[[", 1))
  classes = strsplit(x = unlist(sapply(elements, "[[", 2)), split = "/")
  data$Class = unlist(sapply(classes, "[[", 1))
  data$Cluster = sub(pattern = "_.*", replacement = "", x = data$Subfamily)
  
  return(data)
}
aves = clean.landscape(data = aves)

# get species name
aves$Species = sub(pattern = ".k2p.noCpG.size", replacement = "", x = filename)

# discard entries for satellites not in the list given
boo = aves$Cluster %in% sats
aves = aves[boo,]

data = aves
data$Divergence = data$Divergence * 100
data$RoundDiv = floor(data$Divergence)
### FACTOR IS THE IDENTIFIER YOU CAN WORK ON TO USE GENERAL OR SPECIFIC CLASSES OF REPEATS
data$Factor = paste(data$Cluster, data$RoundDiv, sep = "$")
data_bps = aggregate(Size ~ Factor, data, sum)
data_bps$Cluster = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 1)
data_bps$Divergence = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 2)
size = genomeSize$GenomeSize[genomeSize$Species == unique(data$Species)[1]]
data_bps$Mb = (data_bps$Size / size) * 100
threshold = 0
data_bps_sub = data_bps[data_bps$Mb > threshold,]
data_bps_sub = as.data.frame(data_bps_sub)
data_bps_sub$Species = data$Species[1]
aves_bps = data_bps_sub

aves_bps$Divergence = as.integer(aves_bps$Divergence)
boo = aves_bps$Divergence <= 60
aves_bps = aves_bps[boo,]

# make the satellite column as a factor with all the levels as listed in the sat list
# so that all the landscapes for evry species will have the same legend and colors
# order satellites
levels = as.character(sats)
#levels = as.character(unique(aves_bps$Cluster))
aves_bps$Cluster = factor(x = aves_bps$Cluster, levels = levels[order(nchar(levels), levels)])

# assign a color for each satellite sequence
require(RColorBrewer)
palette = colorRampPalette(brewer.pal(11, "Spectral"))(length(sats))
data_col = data.frame(Satellite = levels[order(nchar(levels), levels)], Colors = palette)
data_col$Colors = as.character(data_col$Colors)
o = match(table = data_col$Satellite, x = levels(droplevels(aves_bps$Cluster)))
colors = data_col$Colors[o]

filename = paste0(sub(pattern = ".k2p.noCpG.size", replacement = "", x = filename),"_divergence_abundance.txt")
write.table(x = aves_bps, file = filename, sep ="\t", quote = F, row.names = F, col.names = T)

# create plot
plot = ggplot(data = aves_bps, aes(x = as.integer(Divergence), y = Mb, fill = Cluster)) + geom_bar(stat = "identity", colour = "gray90", size = 0.15) + xlab("Divergence (%)") + ylab("Proportion of genome (%)") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + xlim(NA, 40) + scale_fill_manual(values = colors) + ggtitle(aves_bps$Species[1])
filename = paste0(plotname, ".pdf")
ggsave(filename = filename, plot = plot, device = "pdf", width = 40, height = 27, units = "cm", scale = .5)
