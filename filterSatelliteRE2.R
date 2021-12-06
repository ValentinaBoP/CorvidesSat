#!/usr/bin/env Rscript
# Usage: Rscript --vanilla filterSatelliteRE2.R satellite.fasta file.blast.out prefix_output

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("Three arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args) == 3) {
  fasta = args[1]
  input = args[2]
  output = args[3]
}

library(data.table)
library(dplyr)

# read input
blast = fread(input, header = FALSE)

# assign column names
names(blast) = c("qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen", "sstrand", "evalue", "length", "pident")

# filter self alignments
boo = blast$qseqid == blast$sseqid
blast = blast[!boo,]

# discard low identities
boo = blast$pident > 85
blast = blast[boo,]

# calculate percentage of aligned sequence
blast$palign = (blast$length / blast$qlen) * 100

# discard short alignments
boo = blast$palign > 80
blast = blast[boo,]

# write blast table filtered
write.table(x = blast, file = paste0(output, ".filtered.blast.out"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# get clusters
clusters = blast %>% group_by(qseqid) %>% summarize(Groups = paste(unique(c(qseqid, sseqid)), collapse = ','))

# write cluster names to file
sink(paste0(output, ".cluster.list"))
print(unique(clusters$Groups))
cat("\n")
sink()

system("mkdir Clusters")
setwd("./Clusters")

for(i in 1:length(unique(clusters$Groups))){

	# get sequence names
	seq_names = paste0(gsub(x = unique(clusters$Groups)[i], pattern = ",", replacement = "\n"), "\n")

	# save sequence names into a list
	cluster_file = paste0("Cluster_", i, ".txt")
	write.table(x = seq_names, file = cluster_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

	# extract fasta sequences from library
	cmd = paste0("perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl ../", fasta, " list ", cluster_file, " > ", paste0("Cluster_", i, ".fasta"))
	system(cmd)

}

# extract sequences that do not belong to the clusters
cmd = paste0("grep '>' ../", fasta, " | cut -c2-")
all_seqs = read.table(text = system(cmd, intern = TRUE), comment.char = "", stringsAsFactors = FALSE, header = FALSE)[,1]

seqs = unique(c(unique(clusters$qseqid), unique(unlist(strsplit(x = clusters$Groups, split = ",")))))
boo = all_seqs %in% seqs
write.table(x = all_seqs[!boo], file = "no_clusters.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
cmd = paste0("perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl ../", fasta, " list no_clusters.txt > no_clusters.fasta")
system(cmd)



