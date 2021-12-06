#!/usr/bin/env Rscript
# Usage: Rscript --vanilla abundance_proportions.R list_sats GenomeSize_summary.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("3 arguments must be supplied", call.=FALSE)
} else if (length(args) == 2) {
  list_sats = args[1]
  sizefile = args[2]
}

# load libraries
library(data.table)
library(dplyr)
library(tidyr)

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

tot = df %>% group_by(Species) %>% summarise(TotalPerc = sum(Proportion))

write.table(x = tot, file = "Satellite_proportion_from_reads.txt", sep = "\t", quote = FALSE, row.names = FALSE)
