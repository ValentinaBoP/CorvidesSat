#!/usr/bin/env Rscript
# Usage: Rscript --vanilla plotFamiliesRPKM.R rpkm_file

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args) == 1) {
  rpkmfile = args[1]
}

#libraries
library(data.table)
library(ggplot2)
library(dplyr)

# read input file
rpkm = fread(rpkmfile, header = TRUE)

# pre-plot processing
rpkm$Family = sub(pattern = "_.*", replacement = "", x = rpkm$Satellite)

# if multiple subfamilies are present, they are compacted and the lowest value only is kept
rpkm = rpkm %>% group_by(Family) %>% filter(rpkm == min(rpkm))

# set levels for plotting
levels = as.character(unique(rpkm$Family))
rpkm$Family = factor(x = rpkm$Family, levels = levels[order(nchar(levels), levels)])

rpkm_fixed = rpkm
rpkm_fixed$rpkm[rpkm_fixed$Family == "sat1"] = 199

plot = ggplot(data = rpkm_fixed, aes(x = Family, y = rpkm)) + geom_bar(stat = "identity", fill = "#087E8B", color = "black") + theme_bw() + ylab("Reads Per Kilobase Million (RPKM)") + xlab("Satellite sequences") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(plot = plot, filename = "RPKM_plot.pdf", units = "cm", width = 20, height = 20, dpi = 300, device = "pdf")
