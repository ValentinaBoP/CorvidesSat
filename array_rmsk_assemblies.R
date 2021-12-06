#!/usr/bin/env Rscript
# Usage: Rscript --vanilla arrays_rmsk.R file.out list_sats

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("3 arguments must be supplied", call.=FALSE)
} else if (length(args) == 2) {
  list_sats = args[2]
  rmskfile = args[1]
}

# load libraries
library(data.table)
library(dplyr)
library(tidyr)

data = fread(rmskfile)

for(i in 1:nrow(data)){

	data$V7[i] = paste(sort(unique(unlist(strsplit(x = data$V4[i], split = ",")))), collapse = ",")
}

data = data[,c(1:4,7,5,6)]

df = data %>% group_by(V7) %>% summarise(Num = n(), Len = max(V6))

df_sub = df[df$Num > 1 & df$Len > 10000,]
df_sub = df_sub[order(df_sub$Num, decreasing = TRUE),]

write.table(x = data, file = "total_arrays_composition_orientation.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

write.table(x = df_sub, file = "array_types_len.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
