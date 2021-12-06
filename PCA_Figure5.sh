# PCA - Figure 5

ml R_packages/3.6.0

cd /home/vpeona/2021SatCorvides/Intermediate/RE2/PresAbs/Hybrids

CODE=/home/vpeona/2021SatCorvides/Code
Rscript --vanilla $CODE/pca_corvus_hybrids.R proportions.txt

cd /home/vpeona/2021SatCorvides/Intermediate/RE2/PresAbs

Rscript --vanilla $CODE/pca_all_species.R proportions.txt

## Figure 5a
#!/usr/bin/env Rscript
# Usage: Rscript --vanilla pca_all_species.R proportions.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("3 arguments must be supplied", call.=FALSE)
} else if (length(args) == 1) {
  propfile = args[1]
}


library(tidyr)
library(data.table)
library(ggbiplot)

data = fread(propfile, header = TRUE)
new = data[,c(1,2,4)]
new_spread = spread(new, Satellite, Proportion)
new_spread = as.data.frame(new_spread)
sat.pca <- prcomp(new_spread[,c(2:length(new_spread))], center = TRUE,scale. = TRUE)

sat.groups = c(rep("BOP", 4), rep("Corvus", 7), rep("BOP", 12))

plot = ggbiplot(sat.pca, ellipse=TRUE, labels=new_spread$Species, groups=sat.groups)

ggsave(plot = plot, filename = "PCA_all_species.pdf", units = "cm", width = 30, height = 20, dpi = 300, device= "pdf")


## Figure 5b
#!/usr/bin/env Rscript
# Usage: Rscript --vanilla pca_corvus_hybrids.R proportions.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("3 arguments must be supplied", call.=FALSE)
} else if (length(args) == 1) {
  propfile = args[1]
}


library(tidyr)
library(data.table)
library(ggbiplot)

data = fread(propfile, header = TRUE)
new = data[,c(1,2,4)]
new_spread = spread(new, Satellite, Proportion)

sat.pca <- prcomp(new_spread[,c(2:23)], center = TRUE,scale. = TRUE)

sat.groups = c(rep("corCon", 3), rep("corCor", 3), rep("Hybrids", 6))

plot = ggbiplot(sat.pca, ellipse=TRUE, labels=new_spread$Species, groups=sat.groups)
ggsave(plot = plot, filename = "PCA_corvus_hybrids_2.2.pdf", units = "cm", width = 30, height = 20, dpi = 300, device= "pdf")
