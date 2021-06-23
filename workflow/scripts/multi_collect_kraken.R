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
print(input_files)
output          <- snakemake@output #  "Collection/multi_{type}.kraken.Rdata"
samples         <- snakemake@params[["samples"]]
output_files    <- snakemake@params[["out_prefix"]] #"Collection/multi_{type}.kraken."


###################################################################################################
## Get data
###################################################################################################

all_krak <- list()
for(i in 1:length(input_files)){
  all_krak[[samples[i]]] <- read.delim(input_files[[i]],stringsAsFactors = F,header=F)[,c(2,4,6)]
  colnames(all_krak[[samples[i]]]) <- c("reads","rank","taxon")
  all_krak[[samples[i]]]$taxon <- sub("^ +","",all_krak[[samples[i]]]$taxon)
}

all_krak <- rbindlist(lapply(1:length(all_krak),function(x) data.frame(all_krak[[x]],
                                                                      "sample"=names(all_krak)[x],
                                                                      stringsAsFactors = F)))

all_krak_ranks <- lapply(unique(all_krak$rank),function(x) all_krak[all_krak$rank==x,c("sample","taxon","reads")])
names(all_krak_ranks) <- unique(all_krak$rank)
all_krak_ranks <- lapply(all_krak_ranks,function(x) tapply(x$reads,list(x$taxon,x$sample),sum))
all_krak_ranks <- lapply(all_krak_ranks,function(x) {
                      x[is.na(x)] <- 0
                      x})

###################################################################################################
## Output
###################################################################################################

for(r in names(all_krak_ranks)){
  write.table(all_krak_ranks[[r]],paste0(output_files[[1]],r,".tsv"),
              col.names=NA,quote = F,sep="\t")
}
save.image(output[[1]])

