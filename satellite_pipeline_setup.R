#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
config = args[1]

        # libraries and functions required for the pipeline
        require(Biostrings)
        require(GenomicRanges)
        source("/home/vpeona/snic2018-8-266/private/Code/SatellitePipeline/satellite_pipeline_functions.R")

	#assign variables
	genome = args[1]
	flank = as.integer(args[2])
	output = args[3]

	# additional variables
	index = paste(genome, ".fai", sep = "")

        # read genome fasta file
        FASTA = readDNAStringSet(filepath = genome, format = "fasta")

        # annotate gaps in the genome
        GAPS = annotateGaps(FASTA = FASTA)
	
		if(nrow(GAPS) != 0){
			
			# get the flanking regions (long as flank) of the gaps
        		range_gaps = getFlankingRangeGaps(GAPS = GAPS, flank = flank)
			range_gaps = as.data.frame(range_gaps)
			range_gaps = range_gaps[,-c(4,5,6)]
		        # transform txt in bed file
        		range_gaps[,1] = as.character(range_gaps[,1])
		        simplifiedScaffolds = sapply(strsplit(x = range_gaps[,1], split = " "), "[[", 1)
		        gap_names = paste(simplifiedScaffolds, "_gaps_", 1:length(simplifiedScaffolds), sep = "")
		        BED = range_gaps[,c(1,2,3)]
   			BED[,4] = gap_names

		} else {
			
			print("There are no gaps in your genome")
			BED = data.frame()
			
		}


        # get the starting and ending regions of scaffolds/contigs
        flanking_contigs = getFlankingRegions(index = index, flank = flank)
		
	BED = rbind(BED, flanking_contigs)

        # fix the length of the start and end coordinates in the bed file (they must not exceed the length of the chromosome
        BED2 = fixBedChrLength(BED, INDEX = index)

        # fix values below 0
        boo = BED2[,2] < 0
        BED2[boo,2] = 0

        # save the new fixed bed file
        write.table(x = BED2, file = "temp_fixed.bed", quote = F, row.names = F, col.names = F, sep = "\t")

        # run bedtools getfasta to retrieve the gap flanking regions
        system(paste("module load bioinfo-tools BEDTools/2.27.1 && bedtools getfasta -fi ", genome, " -bed temp_fixed.bed -fo ", output, " -name", sep = ""))
