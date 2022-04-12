# Usage: Rscript --vanilla genomeSimulation.R genome.fasta genome.fasta.fai sat.lib <perc of sat> <prefix>
# it requires BEDtools and samtools

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 5) {
  stop("Four arguments must be supplied", call.=FALSE)
} else {
  
  fasta = args[1]
	fai = args[2]
	lib = args[3]
	percentage_sat = as.numeric(args[4])
	prefix = args[5]

}

library(dplyr)
library(data.table)

# read index file
index = fread(fai)

# read satellite monomers
cmd = paste0("grep '>' ", lib, " | cut -c2- ")
names = read.table(text = system(cmd, intern = TRUE), stringsAsFactors = FALSE, comment = "")
cmd = paste0("grep -v '>' ", lib)
seqs = read.table(text = system(cmd, intern = TRUE), stringsAsFactors = FALSE)

# data frame with sat names and sat monomers
satellites = data.frame(Repeat = names$V1, Sequence = seqs$V1)
rm(names, seqs)

# set tot percentage of sat to add
genomeSize = sum(index$V2)
bps_sat = floor((genomeSize / 100) * percentage_sat)

# test with 1kb threshold

# counter that needs to be updated each cycle
counter = 1
# I need a while loop until the percentage of sat is reached
# need to set the percentage, therefore the number of bases before starting the insertions

bps_occupied = 0

# name for new set of satellites .fasta
sat_filename = paste0(prefix, "_sat_arrays.fasta")

# source for simulation
source_filename = paste0(prefix, "_genome_simulation_source.fasta")

# name for the simulated output genome fasta file
outfile = paste0(prefix, "_simulated_genome.fasta")

# remove the sat_filename from the folder if any already exists
#system(paste0("rm ", sat_filename))

# BED file with all the insertions to add
master_bed = data.frame(Chromosome = index$V1, Start = 0, End = index$V2)

while(bps_occupied <= bps_sat){
	# get random chromosome
	chr = sample(index$V1, 1)
 
	# get random location on chr
	temp = master_bed[master_bed$Chromosome == chr,][sample(1:nrow(master_bed[master_bed$Chromosome == chr,]), 1),]
	new_locus = sample(temp[,2]:temp[,3], 1)
	rm(temp)
	#new_locus = sample(master_bed$Start[master_bed$Chromosome == chr]:master_bed$End[master_bed$Chromosome == chr], 1)

	# create sequence to insert in the new locus
	sat = sample(1:nrow(satellites), 1)
	n_monomers = sample(1:100, 1)
	new_array = paste(rep(satellites$Sequence[sat], n_monomers), collapse = "")
	new_array_name = paste(chr, satellites$Repeat[sat], n_monomers, counter, sep = "_")

	# divide the chr where the new_locus is
	# locate the right interval to cut
	row = which(master_bed$Chromosome == chr & master_bed$Start <= new_locus & master_bed$End >= new_locus)
	new_5end = data.frame(Chromosome = chr, Start = master_bed$Start[row], End = new_locus)
	new_3end = data.frame(Chromosome = chr, Start = new_locus, End = master_bed$End[row])
	new_locus_line = data.frame(Chromosome = new_array_name, Start = 0, End = nchar(new_array))

	# add the rows to the data.frame
	master_bed = add_row(master_bed[-row,], .after = row-1, new_5end)
	master_bed = add_row(master_bed, .after = row, new_locus_line)
	master_bed = add_row(master_bed, .after = row + 1, new_3end)

	counter = counter + 1

	bps_occupied = bps_occupied + nchar(new_array)

	sink(sat_filename, append = TRUE)
	cat(paste0(">", new_array_name, "\n", new_array, "\n"))
	sink()

}

print(paste("bps occupied: ", bps_occupied))

# create the new genome one chromosome at the time

# concatenate the initial genome with the new satellite arrays
cmd = paste0("cat ", fasta, " ", sat_filename, " > ", source_filename)
system(cmd)

for(i in 1:nrow(index)){

	# get chromosome name and subset master_bed for that chromosome
	chr = index$V1[i]
	chr_bed = master_bed[grepl(pattern = paste0(chr, "$|", chr, "_"), x = master_bed$Chromosome),]

	# save the new BED file
	bed_filename = paste0(chr, ".bed")
	write.table(x = chr_bed, file = bed_filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

	# get sequences for the new chr from genome_simulation_source.fasta
	fasta_filename = paste0(chr, "2stitch.fasta")
	cmd = paste("bedtools getfasta -fi", source_filename, "-bed", bed_filename, "-fo", fasta_filename)
	system(cmd)

	cmd = paste0("grep -v '>' ", fasta_filename, " | tr -d '\\n' | sed '1s/^/>", chr,"\\n/' | sed -e '$a\\' >> ", outfile)
	system(cmd)

	cmd = paste("rm", bed_filename, fasta_filename)
	system(cmd)

}





