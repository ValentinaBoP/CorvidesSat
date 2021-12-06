#!/usr/bin/env Rscript
# Usage: Rscript --vanilla presence_sat.R list_sats GenomeSize_summary.txt figure_name

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("3 arguments must be supplied", call.=FALSE)
} else if (length(args) == 3) {
  list_sats = args[1]
  sizefile = args[2]
  plotname = args[3]
}

# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# read list of the satellites present in all species taken into account
sats = fread(list_sats, header = FALSE)$V1

# read library sizes for each species
sizes = fread(sizefile, header = TRUE)

# get names of the files with sat abundance and names
filenames = list.files(pattern = "abundance.clusters.filter.txt")

# create dataframe that will contain T/F for presence absence of sats listes in list_sats
df = data.frame(Species = character(), Satellite = character(), Value = integer(), Proportion = numeric())

# add presence absence for each sat (0 or 1) per species
for(i in 1:length(filenames)){

	species = sub(pattern = ".abundance.clusters.filter.txt", replacement = "", x = filenames[i])
	data = fread(filenames[i], header = TRUE)
	o = match(x = data$Cluster, table = sats)
	prop = rep(0, length(sats))
	prop[o] = (data$Count / sizes$GenomeSize[sizes$Species == species]) * 100
	newData = data.frame(Species = species, Satellite = sats, Value = as.integer(sats %in% data$Cluster), Proportion = prop)

	df = rbind(df, newData)

}

# pre-plot
# order satellites
levels = as.character(unique(df$Satellite))
df$Satellite = factor(x = df$Satellite, levels = levels[order(nchar(levels), levels)])

# order species
levels = rev(c("corCor", "corCon", "corBra", "corWoo", "corMon", "corSpl", "corDau", "cinMag", "cinReg", "parRag", "parRub", "astSte", "astRot", "parBre", "dreAlb", "epiMey", "ptiInt", "ptiMag", "parHel", "parLaw", "manKer", "manCha", "lycPyr"))
df$Species = factor(x = df$Species, levels = levels)

# replace values higher than 1 with 0, they will be colored with Inkscape to help the visualisation of the proportions otherwise flattened by the high values
df$PropFixed = df$Proportion
df[df$PropFixed >= 1,]
df$PropFixed[df$PropFixed >= 1] = 0

write.table(x = df, file = "proportions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# plot
plot = ggplot(data = df, aes(x = Satellite, y = Species)) + geom_tile(aes(fill = Value)) + scale_fill_distiller(palette = "Spectral", direction = -1) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
plotname = paste0(plotname, "_noProp.pdf")
ggsave(plot = plot, filename = plotname, units = "cm", width = 21, height = 30, dpi = 300, device = "pdf")

# plot heatmap by proportion of genome
#plot = ggplot(data = df, aes(x = Satellite, y = Species)) + geom_tile(aes(fill = Proportion)) + scale_fill_distiller(palette = "OrRd", direction = -1) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
plot = ggplot(data = df, aes(x = Satellite, y = Species)) + geom_tile(aes(fill = PropFixed)) + scale_fill_gradient(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
plotname = paste0(plotname, "_withProp.pdf")
ggsave(plot = plot, filename = plotname, units = "cm", width = 21, height = 30, dpi = 300, device = "pdf")
