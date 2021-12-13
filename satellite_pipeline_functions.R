annotateGaps = function(filepath = "0", FASTA){

  library(Biostrings)

  if(filepath == "0" & class(FASTA) == "DNAStringSet"){

    print("Analyzing DNAStringSet")

  } else if(filepath == "0" & class(FASTA) != "DNAStringSet") {

    print("DNAStringSet object is needed!")
    break()

  } else {

    FASTA = readDNAStringSet(filepath = filepath, format = "fasta")

  }

  GAPS = data.frame(start = integer(), end = integer(), width = integer())

  for(i in 1:length(FASTA)){

    y = maskMotif(FASTA[[i]],'N')
    z = as(gaps(y),"Views")
    NEW = as.data.frame(ranges(z))

    if(nrow(NEW) > 0){

      NEW$scaffold = names(FASTA[i])
      GAPS = rbind(GAPS, NEW)

    } else {

      next()

    }

  }

  return(GAPS)

}

getFlankingRangeGaps = function(GAPS, flank = 500){

  library(GenomicRanges)

  out_gaps = data.frame(start = integer(), end = integer(), scaffold = character(), width = integer())

  NEW = data.frame(start = GAPS$start - flank, end = GAPS$start, scaffold = GAPS$scaffold, width = GAPS$width)
  out_gaps = rbind(out_gaps, NEW)

  NEW = data.frame(start = GAPS$end, end = GAPS$end + flank, scaffold = GAPS$scaffold, width = GAPS$width)
  out_gaps = rbind(out_gaps, NEW)
  out_gaps = out_gaps[order(out_gaps[,3], out_gaps[,1]),]

  range_gaps = GRanges(out_gaps$scaffold, IRanges(out_gaps$start, out_gaps$end), mcols = data.frame(width = out_gaps$width))

  return(range_gaps)
}

getFlankingRegions = function(index, flank){

        index = read.table(file = index, sep = "\t", header = F, stringsAsFactors = F)
        flanking_5 = data.frame(seqnames = index[,1], start = 1, end = flank, V4 = paste(index[,1], "_flank5_", 1:nrow(index), sep = ""))
        flanking_3 = data.frame(seqnames = index[,1], start = (index[,2] - flank), end = index[,2], V4 = paste(index[,1], "_flank3_", 1:nrow(index), sep = ""))
        return(rbind(flanking_5, flanking_3))

}

filterRM = function(RM, gaps, method = "aggressive", out){

                DISCARD = c("Satellite", "Simple_repeat", "Low_complexity")
                RM = RM[!(RM$Family %in% DISCARD),]

                gaps = gaps[gaps[,4] %in% RM$Sequence,]

                write.table(x = gaps, file = out, quote = F, sep = '\t', row.names = FALSE, col.names = FALSE)
        return(gaps)
}
readRM = function(filepath){

  RM = read.table(file = filepath, skip = 3, fill = TRUE, stringsAsFactors = FALSE, row.names = NULL)
  RM = RM[!is.na(RM[,2]),]
  names(RM) = c("Score", "Divergence", "Deletion", "Insertion", "Sequence", "Begin", "End", "SequenceLeft", "Strand", "Repeat", "Family", "RepeatBegin", "RepeatEnd", "RepeatLeft", "ID")
  RM$Strand = sub(pattern = "C", replacement = '-', x = RM$Strand)
  return(RM)

}

parseConfigFile = function(config){

        list_vars = read.table(file = config, header = F, stringsAsFactors = FALSE, sep = "\t")
        return(list_vars)
}

writeRepeatMaskerJob = function(out, flank, dir, species, RMJob, config){

        sink(RMJob)
        cat("#!/bin/bash -l")
        cat("\n")
        cat(paste("#SBATCH -J RM_", out, "_gapsFlanking_", flank, sep = ""))
        cat("\n")
        cat(paste("#SBATCH -o RM_", out, "_gapsFlanking_", flank, ".output", sep = ""))
        cat("\n")
        cat(paste("#SBATCH -e RM_", out, "_gapsFlanking_", flank, ".error", sep = ""))
        cat("\n")
        cat("#SBATCH --mail-user valentina.peona90@gmail.com")
        cat("\n")
        cat("#SBATCH --mail-type=ALL")
        cat("\n")
        cat("#SBATCH -t 01:00:00")
        cat("\n")
        cat("#SBATCH -A snic2019-8-79")
        cat("\n")
        cat("#SBATCH -p core")
        cat("\n")
        cat("#SBATCH -n 8")
        cat("\n")
        cat("cd", dir)
        cat("\n")
        cat(paste("cp ", out, "_gapsFlanking_", flank, "_fixed.fasta $SNIC_TMP", sep = ""))
        cat("\n")
        cat("cd $SNIC_TMP")

        cat("\n")
        cat("module load bioinfo-tools RepeatMasker/4.0.7")
        cat("\n")

        cat(paste("REF=", paste(out, "_gapsFlanking_", flank, "_fixed.fasta", sep = ""), sep = ""))
        cat("\n")

        cat("\n")
        cat(paste("RepeatMasker -pa 8 -a -xsmall -gccalc -species ", species, " $REF", sep = ''))
        cat("\n")

        cat("cp", paste(out, "_gapsFlanking_", flank, "_fixed.fasta.*", sep = ""), dir)
        cat("\n")

        cat(paste("cd", dir))
        cat("\n")

        cat(paste("Rscript --vanilla /crex/proj/snic2018-8-266/private/Code/SatellitePipeline/satellite_pipeline_step2.R ", config, sep = ""))
        cat("\n")

        sink()

}

writeRepeatModelerJob = function(out, flank, dir, database, RModJob, config){

sink(RModJob)
cat("#!/bin/bash -l")
cat("\n")
cat(paste("#SBATCH -J RModeler_", out, "_gapsFlanking_", flank, sep = ""))
cat("\n")
cat(paste("#SBATCH -o RModeler_", out, "_gapsFlanking_", flank, ".output", sep = ""))
cat("\n")
cat(paste("#SBATCH -e RModeler_", out, "_gapsFlanking_", flank, ".error", sep = ""))
cat("\n")
cat("#SBATCH --mail-user valentina.peona90@gmail.com")
cat("\n")
cat("#SBATCH --mail-type=ALL")
cat("\n")
cat("#SBATCH -t 24:00:00")
cat("\n")
cat("#SBATCH -A snic2019-8-79")
cat("\n")
cat("#SBATCH -p core")
cat("\n")
cat("#SBATCH -n 16")
cat("\n")
cat("module load bioinfo-tools RepeatModeler/1.0.8_RM4.0.7")
cat("\n")
cat("cd", dir)
cat("\n")
cat(paste("BuildDatabase -name ", database, " -engine ncbi ", out, "_gapsFlanking_", flank, "_fixed_filter.fasta", sep = ""))
cat("\n")
cat(paste("cp ", database, "* $SNIC_TMP", sep = ""))
cat("\n")
cat("cd $SNIC_TMP")
cat("\n")
cat(paste("RepeatModeler -pa 16 -database", database))
cat("\n")
cat(paste("mkdir ", dir, "/RModeler", sep = ""))
cat("\n")
cat("mv $SNIC_TMP/RM_*/consensi.*", paste(dir, "/RModeler", sep = ""))
cat("\n")

cat(paste("cd ", dir, sep = ""))
cat("\n")

cat(paste("Rscript --vanilla /crex/proj/snic2018-8-266/private/Code/SatellitePipeline/satellite_pipeline_step3.R ", config, sep = ""))
cat("\n")

sink()

}

writeRPipeJob = function(out, flank, dir, home_dir, RPipeJob){

sink(RPipeJob)
cat("#!/bin/bash -l")
cat("\n")
cat(paste("#SBATCH -J rmodelerPipeline_", out, "_gapsFlanking_", flank, sep = ""))
cat("\n")
cat(paste("#SBATCH -o rmodelerPipeline_", out, "_gapsFlanking_", flank, ".output", sep = ""))
cat("\n")
cat(paste("#SBATCH -e rmodelerPipeline_", out, "_gapsFlanking_", flank, ".error", sep = ""))
cat("\n")
cat("#SBATCH --mail-user valentina.peona90@gmail.com")
cat("\n")
cat("#SBATCH --mail-type=ALL")
cat("\n")
cat("#SBATCH -t 24:00:00")
cat("\n")
cat("#SBATCH -A snic2019-8-79")
cat("\n")
cat("#SBATCH -p core")
cat("\n")
cat("#SBATCH -n 16")
cat("\n")
cat("module load bioinfo-tools blast/2.2.31+")
cat("\n")
cat("module load bioinfo-tools MAFFT/7.310")
cat("\n")
cat("unset $MAFFT_BINARIES")
cat("\n")

cat(paste("CONS=", dir, "/RModeler/consensi.fa.classified", sep = ""))
cat("\n")
cat(paste("REF=", genome, sep = ""))
cat("\n")
cat(paste("DIR=", dir, "/RModeler", sep = ""))
cat("\n")
cat("cd $DIR")
cat("\n")

cat("makeblastdb -in $REF -out $REF -dbtype nucl -parse_seqids")
cat("\n")

cat("awk '{print$1;}' $CONS > consensi.fa.classified.clean")
cat("\n")

cat("CONS=consensi.fa.classified.clean")
cat("\n")

cat("sed -i 's/\\//_/g' $CONS")
cat("\n")

cat("perl /crex/proj/snic2018-8-266/private/Code/SatellitePipeline/repeatModelerPipeline4.pl $REF $REF $CONS")
cat("\n")

cat("module purge")
cat("\n")
cat("module load bioinfo-tools T-Coffee/11.00.8cbe486")
cat("\n")
cat(paste("export CACHE_4_TCOFFEE=", home_dir, sep = ""))
cat("\n")

cat("rm *emp.out")
cat("\n")
cat("mkdir final")
cat("\n")

cat("cd aligned; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d \".\"`; cat $i | perl -ne 'chomp;s/>\\s+/>/;if(/>(\\S+)/){$id{$1}++;$id2=$1;}if($id{$id2}==1){print \"$_\\n\"}' >../final/$name.fa; done; cd ../")
cat("\n")

cat("cd final; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d \".\"`; t_coffee -other_pg seq_reformat -in $i -action +rm_gap 95 >$name.gaps95.fa; done; cd ../")
cat("\n")
sink()

}

fixBedChrLength = function(BED, INDEX){

INDEX = read.table(INDEX, stringsAsFactors=F, header = F, sep = "\t")

BED[,1] = as.character(BED[,1])
BED[,1] = sapply(strsplit(x = BED[,1], split = " "), "[[", 1)

ELEMENT = INDEX[,1]
for(i in 1:length(ELEMENT)){

        BED[BED[,1] == ELEMENT[i] & BED[,3] > INDEX[i,2],3] = INDEX[i,2]

}

return(BED)
}
