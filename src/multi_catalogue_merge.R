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
library(data.table)

###################################################################################################
## Files from snakemake
###################################################################################################

input_files     <- snakemake@input #actual files,determined based on step and what's there
output          <- snakemake@output #  "Catalogue/catalogue.tsv"
krak            <- snakemake@config[["catalogue"]][["kraken"]][["do"]]
map             <- snakemake@config[['catalogue']][['mapping']][['do']]
types           <- snakemake@params[["types"]]

###################################################################################################
## Get data
###################################################################################################

cat <- read.delim(input_files[[1]],stringsAsFactors = F)
cat <- cat[,c(1,ncol(cat))]
colnames(cat)[2] <- "members" 
cat <- merge(cat,read.delim(input_files[[2]],stringsAsFactors = F),by="gene",all=T)

if(krak){
  n <- 4
  krakenTab <- read.delim(input_files[[3]],stringsAsFactors = F)
  krakenTab$tax <- apply(krakenTab,1,function(x) paste(x[!is.na(x)],sep=";",collapse=";"))
  cat <- merge(cat,krakenTab,by.x="gene",by.y="contig",all=T)
  rm(krakenTab)
}else{
  n <- 3
}

if(map){
  for(type in types){
    fcfile <- read.delim(input_files[[n]],comment.char = "#",stringsAsFactors = F)
    fcfile <- fcfile[,c(1,7:ncol(fcfile))]
    cat <- merge(cat,fcfile, by.x="gene",by.y="Geneid",all=T)
    rm(fcfile)
    n <- n + 1
  }
}

###################################################################################################
## Output
###################################################################################################

write.table(cat,output[[1]],
              col.names=NA,quote = F,sep="\t")

