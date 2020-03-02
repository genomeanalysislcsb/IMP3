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

input_files         <- snakemake@input #actual files,determined based on step and what's there
print(input_files)
output             <- snakemake@output # "Collection/multi_{type}.{db}.tsv",category="Collection"),
                                        # "Collection/multi_{type}.{db}.RDS"
samples <- snakemake@params[["samples"]]

###################################################################################################
## Get data
###################################################################################################

all_feat <- list()
for(i in 1:length(input_files)){
  all_feat[[samples[i]]] <- read.delim(input_files[[i]],comment.char = "#",stringsAsFactors = F)[,c(1,7)]
  colnames(all_feat[[samples[i]]])[ncol(all_feat[[samples[i]]])] <- "reads"
}

all_feat <- rbindlist(lapply(1:length(all_feat),function(x) data.frame(all_feat[[x]],
                                                                      "sample"=names(all_feat)[x],
                                                                      stringsAsFactors = F)))
all_feat <- tapply(all_feat$reads,list(all_feat$Geneid,all_feat$sample),sum)
all_feat[is.na(all_feat)] <- 0

###################################################################################################
## Output
###################################################################################################

write.table(all_feat,output[[1]],col.names=NA,quote = F,sep="\t")
saveRDS(all_feat,output[[2]])

