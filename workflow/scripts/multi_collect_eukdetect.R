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
output          <- snakemake@output[[1]] #  "Collection/multi_{type}.eukdetect.Rdata"
samples         <- snakemake@params[["samples"]]
print(samples)
output_files    <- snakemake@output[-1] #expand("Collection/multi_{type}.eukdetect.{kind}.tsv",kind=["Read_counts","Total_marker_coverage","Percent_identity"])

###################################################################################################
## Get data
###################################################################################################

all_euk <- list()
for(i in 1:length(input_files)){
  all_euk[[samples[i]]] <- read.delim(input_files[[i]],stringsAsFactors = F)
  all_euk[[samples[i]]]$Total_marker_coverage <- as.numeric(gsub("%$","",all_euk[[samples[i]]]$Total_marker_coverage))
  all_euk[[samples[i]]]$Percent_identity <- as.numeric(gsub("%$","",all_euk[[samples[i]]]$Percent_identity))
 # colnames(all_euk[[samples[i]]]) <- c("taxstrin","reads")
 # all_euk[[samples[i]]]$rank <- sapply(all_mOTU[[samples[i]]]$taxstring, 
 #                                       function(x) {
 #                                         tax_vec <- unlist(strsplit(x,split="|",fixed=T))
 #                                         gsub("__.+","",tax_vec[length(tax_vec)])})
}

all_euk <- rbindlist(lapply(1:length(all_euk),function(x) data.frame(all_euk[[x]],
                                                                      "sample"=names(all_euk)[x],
                                                                      stringsAsFactors = F)))

#all_euk_kinds <- lapply(c("Read_counts","Total_marker_coverage","Percent_identity"),
#                         function(x) all_euk[,c("Name",x)])

all_euk_kinds <- lapply(c("Read_counts","Total_marker_coverage","Percent_identity"),
                         function(x) tapply(all_euk[,get(x)],list(all_euk$Name,all_euk$sample),sum))
all_euk_kinds <- lapply(all_euk_kinds,function(x) {
                        x[is.na(x)] <- 0
                        x})

###################################################################################################
## Output
###################################################################################################
for(i in 1:length(all_euk_kinds)){
  r <- c("Read_counts","Total_marker_coverage","Percent_identity")[i]
  out <- grep(paste0("eukdetect.",r,".tsv"),unlist(output_files),value=T)
  write.table(all_euk_kinds[[i]],out,
              col.names=NA,quote = F,sep="\t")
}
save.image(output[[1]])

