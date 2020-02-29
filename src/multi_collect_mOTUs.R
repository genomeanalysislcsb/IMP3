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
output          <- snakemake@output[[1]] #  "Collection/multi_{type}.mOTUs.Rdata"
samples         <- snakemake@params[["samples"]]
output_files    <- snakemake@output[-1] #expand("Collection/multi_{type}.mOTUs.{rank}.tsv",rank=["k","p","c","o","f","g","s"])

###################################################################################################
## Get data
###################################################################################################

all_mOTU <- list()
for(i in 1:length(input_files)){
  all_mOTU[[samples[i]]] <- read.delim(input_files[i],stringsAsFactors = F)
  colnames(all_mOTU[[samples[i]]]) <- c("taxstring","reads")
  all_mOTU[[samples[i]]]$rank <- sapply(all_mOTU[[samples[i]]]$taxstring, 
                                        function(x) {
                                          tax_vec <- unlist(strsplit(x,split="|",fixed=T))
                                          gsub("__.+","",taxvec[length(taxvec)])})
}

all_mOTU <- rbindlist(lapply(1:length(all_mOTU),function(x) data.frame(all_mOTU[[x]],
                                                                      "sample"=names(all_mOTU)[x],
                                                                      stringsAsFactors = F)))

all_mOTU_ranks <- lapply(c("k","p","c","o","f","g","s"),
                         function(x) all_mOTU[all_mOTU$rank==x,c("sample","taxstring","reads")])
all_mOTU_ranks <- lapply(all_mOTU_ranks,function(x) tapply(x,list(x$reads,x$sample),sum))
all_mOTU_ranks <- lapply(all_mOTU_ranks,function(x) x[is.na(x)] <- 0)

###################################################################################################
## Output
###################################################################################################

for(r in names(all_mOTU_ranks)){
  out <- grep(paste0("mOTUs.",r,".tsv"),value=T)
  write.table(all_mOTU_ranks[[r]],out,
              col.names=NA,quote = F,sep="\t")
}
save.image(output[[1]])

