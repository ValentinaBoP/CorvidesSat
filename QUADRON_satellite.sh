cd /home/vpeona/2021SatCorvides/Intermediate
mkdir Quadron
cd Quadron

ln -s ../../Results/corvides_sat1.4.fasta .
ml R_packages/3.5.0

CODE=/home/vpeona/2021SatCorvides/Code

Rscript --vanilla $CODE/quadronPerChromosome.R corvides_sat1.4.fasta




#### R

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(seqinr)
load("/home/vpeona/Quadron/Quadron.lib")

fasta = args[1]

seqs = read.fasta(fasta, seqonly = TRUE)
names = unlist(attributes(read.fasta(fasta)))

for(i in 1:length(names)){
	
	filename = paste(names[i], ".fasta", sep = "")
	write.fasta(sequences = seqs[i], names = names[i], file.out = filename)
	
	out = paste(filename, ".Quadron", sep = "")
	Quadron(FastaFile = filename, OutFile = out, nCPU = 1)

}
