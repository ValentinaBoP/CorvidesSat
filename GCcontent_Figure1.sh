## GC content calculation for consensus sequences - Figure 1
ml R_packages/3.6.0

WD=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Intermediate/GCcontent
mkdir -p $WD && cd $WD

CODE=/proj/sllstore2017073/private/Valentina/2021SatCorvides/Code

ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.3.fasta .

Rscript --vanilla $CODE/gc_content.R corvides_sat1.3.fasta Figure_gc_content_1.0

Rscript --vanilla $CODE/gc_content.R gallus.fasta Figure_gc_content_gallus_1.0

ln -s /proj/sllstore2017073/private/Valentina/2021SatCorvides/Results/corvides_sat1.4.fasta .

Rscript --vanilla $CODE/gc_content.R corvides_sat1.4.fasta Figure_gc_content_2.0

###
#!/usr/bin/env Rscript
# Usage: Rscript --vanilla gc_content.R seq.fasta prefix

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args) == 2) {
  filename = args[1]
  prefix = args[2]
}

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(seqinr)

# input file sequences
fasta = unlist(read.fasta(filename, seqonly = TRUE))

# input sequence names
fasta_names = unlist(attributes(read.fasta(filename)))

# data frame with results
df = data.frame(Sequence = character(), GC = numeric(), AT = numeric())

for(i in 1:length(fasta)){

  # convert in vector of letters
  seq = s2c(fasta[i])

  # calculate exact GC and AT content
  gc = GC(seq, exact = TRUE)
  at = 1 - gc

  # put results in the data frame
  newLine = data.frame(Sequence = fasta_names[i], GC = gc, AT = at)
  df = rbind(df, newLine)

}

# proportion to percentage
df$GC = df$GC * 100
df$AT = df$AT * 100

# assign BOP or crow
boo = grepl(pattern = "cor|crow", x = df$Sequence)
df$Group[boo] = "Corvus"
df$Group[!boo] = "BOP"
df$Group = factor(x = df$Group, levels = c("BOP", "Corvus"))

plot = ggplot(data = df, aes(x = GC, fill = Group)) + geom_histogram(binwidth = 3, color = "black", alpha = 0.5) + theme_bw() + xlab("G+C content (%)") + ylab("Number of satellite DNA consensus sequences") + scale_fill_manual(values = c("#087E8B", "#F2BB05")) + geom_vline(xintercept = 43, linetype = "dotted", color = "black", size = 1)
ggsave(plot = plot, filename = paste0(prefix, ".pdf"), width = 20, height = 20, dpi = 300, device = "pdf", units = "cm")

write.table(x = df, file = "GC_content_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

