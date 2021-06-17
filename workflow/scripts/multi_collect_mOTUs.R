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
output          <- snakemake@output[[1]] #  "Collection/multi_{type}.mOTUs.Rdata"
samples         <- snakemake@params[["samples"]]
print(samples)
output_files    <- snakemake@output[-1] #expand("Collection/multi_{type}.mOTUs.{rank}.tsv",rank=["k","p","c","o","f","g","s"])

###################################################################################################
## Get data
###################################################################################################

all_mOTU <- list()
for(i in 1:length(input_files)){
  all_mOTU[[samples[i]]] <- read.delim(input_files[[i]],stringsAsFactors = F)
  colnames(all_mOTU[[samples[i]]]) <- c("taxstring","reads")
  all_mOTU[[samples[i]]]$rank <- sapply(all_mOTU[[samples[i]]]$taxstring, 
                                        function(x) {
                                          tax_vec <- unlist(strsplit(x,split="|",fixed=T))
                                          gsub("__.+","",tax_vec[length(tax_vec)])})
}

all_mOTU <- rbindlist(lapply(1:length(all_mOTU),function(x) if(nrow(all_mOTU[[x]])>0) data.frame(all_mOTU[[x]],
                                                                      "sample"=names(all_mOTU)[x],
                                                                      stringsAsFactors = F)))

all_mOTU_ranks <- lapply(c("k","p","c","o","f","g","s"),
                         function(x) all_mOTU[all_mOTU$rank==x,c("sample","taxstring","reads")])
all_mOTU_ranks <- lapply(all_mOTU_ranks,function(x) tapply(x$reads,list(x$taxstring,x$sample),sum))
all_mOTU_ranks <- lapply(all_mOTU_ranks,function(x) {
                        x[is.na(x)] <- 0
                        x})

###################################################################################################
## Output
###################################################################################################
for(i in 1:length(all_mOTU_ranks)){
  r <- c("k","p","c","o","f","g","s")[i]
  out <- grep(paste0("mOTUs.",r,".tsv"),unlist(output_files),value=T)
  write.table(all_mOTU_ranks[[i]],out,
              col.names=NA,quote = F,sep="\t")
}
save.image(output[[1]])

