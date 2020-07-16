#! /usr/bin/prun R-3.4.0-0 Rscript

args <- commandArgs(trailingOnly=TRUE)

##### load packages
library(methods)
library(stats)
library(edgeR)
library(plyr)
library(data.table)

### load input files
outpath <- args[1]
samplename <- args[2]
reads <- as.character(args[3])
out <- as.character(args[4])


### Keep only not misleading file
data <- fread(file = reads, sep = "\t")
sum <- ddply(data, .variables = .(V14, V12, V20), transform, 
             misleading = length(unique(V23))>1)
misleading <- sum[sum$misleading == TRUE,]
remove_misleading <- write.table(misleading$V1, sep = "\n", file = out, 
                                 quote = FALSE, row.names = FALSE, col.names = FALSE)