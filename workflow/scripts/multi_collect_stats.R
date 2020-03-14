#!/bin/R
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

###################################################################################################
## Load required packages
###################################################################################################
libs <- paste0(Sys.getenv("CONDA_PREFIX"),"/lib/R/library")
.libPaths(libs)
.libPaths()
library(RColorBrewer)
library(data.table)

###################################################################################################
## Files from snakemake
###################################################################################################

input_files     <- snakemake@input #actual files,determined based on step and what's there
output          <- snakemake@output #  "Collection/multi_stats.tsv","Collection/multi_stats.RDS"
samples         <- snakemake@params[["samples"]]

###################################################################################################
## Get data
###################################################################################################

all_stats <- list()
for(i in 1:length(input_files)){
  all_stats[[samples[i]]] <- read.delim(input_files[[i]],stringsAsFactors = F)
  all_stats[[samples[i]]]$feature <- apply(all_stats[[samples[i]]][,1:2],1,
                                           function(x) paste(x,sep=":",collapse=":"))
}

all_stats <- rbindlist(lapply(1:length(all_stats),function(x) data.frame(all_stats[[x]][,3:4],
                                                                      "sample"=names(all_stats)[x],
                                                                      stringsAsFactors = F)))

all_stats <- tapply(all_stats$number,list(all_stats$feature,all_stats$sample),sum)


###################################################################################################
## Output
###################################################################################################

write.table(all_stats,output[[1]],
              col.names=NA,quote = F,sep="\t")
saveRDS(all_stats,output[[2]])

