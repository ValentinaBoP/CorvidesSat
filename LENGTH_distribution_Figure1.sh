#!/bin/bash -l 
#SBATCH -J LENGTH_distribution
#SBATCH -o LENGTH_distribution.output
#SBATCH -e LENGTH_distribution.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 00:15:00
#SBATCH -A snic2020-5-680
#SBATCH -p core
#SBATCH -n 1


# Figure length distribution
ml bioinfo-tools samtools/1.12 R_packages/3.6.0

CODE=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Code

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/Distribution
mkdir -p $WD && cd $WD

ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4.fasta .

samtools faidx corvides_sat1.4.fasta

# divide the sequences in lengths shorter than 1000 and longer than 1000 bps for plotting reasons
awk '$2 <= 1000 {print $1, $2}' corvides_sat1.4.fasta.fai > corvides_sat1.4_shorter1000.fai
awk '$2 >= 1000 {print $1, $2}' corvides_sat1.4.fasta.fai > corvides_sat1.4_longer1000.fai
awk '{if($1 ~ /crowSat1#Satellite/) {print $1, 1400} else {print $0}}' OFS="\t" corvides_sat1.4_longer1000.fai > temp
mv temp corvides_sat1.4_longer1000.fai

Rscript --vanilla $CODE/length_distribution_barplot.R corvides_sat1.4_shorter1000.fai 50 Figure_lenght_shorter_1.4
Rscript --vanilla $CODE/length_distribution_barplot.R corvides_sat1.4_longer1000.fai 500 Figure_lenght_longer_1.4

###
#!/usr/bin/env Rscript
# Usage: Rscript --vanilla length_distribution_barplot.R index binwidth prefix

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args) == 3) {
  indexname = args[1]
  prefix = args[3]
  binwidth = as.integer(args[2])
}

# libraries
library(ggplot2)
library(data.table)

# read input file
index = fread(indexname, header = FALSE)
names(index) = c("Consensus", "Length")

# assign BOP or crow
boo = grepl(pattern = "cor|crow", x = index$Consensus)
index$Group[boo] = "Corvus"
index$Group[!boo] = "BOP"
index$Group = factor(x = index$Group, levels = c("BOP", "Corvus"))

plot = ggplot(data = index, aes(x = Length, fill = Group)) + geom_histogram(binwidth = binwidth, color = "black", alpha = 0.5) + xlab("Monomer size (bps)") + ylab("Number of satDNA consensus sequences") + theme_bw() + ylim(NA,38) + scale_fill_manual(values = c("#087E8B", "#F2BB05"))
ggsave(plot = plot, filename = paste0(prefix, ".pdf"), width = 20, height = 20, dpi = 300, device = "pdf", units = "cm")
