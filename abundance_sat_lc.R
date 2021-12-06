#!/usr/bin/env Rscript
# Usage: Rscript --vanilla abundance_sat_lc.R file.k2p.noCpG

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args) == 1) {
  filename = args[1]
}

library(data.table)
library(dplyr)

# read input
data = fread(filename)

# read index file of the satellite library
index = fread("../../../../Results/corvides_sat1.3_dimers.fasta.fai")

# filter for satellite
boo = grepl(pattern = "Satellite", x = data$V10)
data = data[boo,]

# calculate abundance for each satellite
data_count = data %>% group_by(V10) %>% summarise(Count = sum(V17)) %>% arrange(desc(Count))

# assign lengths
o = match(table= index$V1, x = data_count$V10)
data_count$Length = index$V2[o]/2

# write table
fileout = sub(pattern = "k2p.noCpG.size", replacement = "abundance.txt", x = filename)
write.table(x = data_count, file = fileout, sep = "\t", quote = F, row.names = F, col.names = TRUE)

# write table with filter for abundance
boo = data_count$Count > 50000
data_count = data_count[boo,]
fileout = sub(pattern = "k2p.noCpG.size", replacement = "abundance.filter.txt", x = filename)
write.table(x = data_count, file = fileout, sep = "\t", quote = F, row.names = F, col.names = TRUE)

# aggregate by cluster and calculate abundance
data$Cluster = sub(pattern = "_.*", replacement = "", x = data$V10)
o = match(table= index$V1, x = data$V10)
data$Length = index$V2[o]/2
data_count_cluster = data %>% group_by(Cluster) %>% summarise(Count = sum(V17), Length = paste(unique(Length), collapse = ",")) %>% arrange(desc(Count))

# write table
fileout = sub(pattern = "k2p.noCpG.size", replacement = "abundance.clusters.txt", x = filename)
write.table(x = data_count_cluster, file = fileout, sep = "\t", quote = F, row.names = F, col.names = TRUE)

# write table with filter for abundance
boo = data_count_cluster$Count > 50000
data_count_cluster = data_count_cluster[boo,]
fileout = sub(pattern = "k2p.noCpG.size", replacement = "abundance.clusters.filter.txt", x = filename)
write.table(x = data_count_cluster, file = fileout, sep = "\t", quote = F, row.names = F, col.names = TRUE)
