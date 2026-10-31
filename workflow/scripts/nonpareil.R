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
library(Nonpareil)

###################################################################################################
## Files from snakemake
###################################################################################################

input_data          <- snakemake@input[[1]] #"Stats/mg/mg.{read}_nonpareil.npo"
outplot             <- snakemake@output[[1]] #Visualization/mg.{read}_nonpareil.curve.png"
outtab              <- snakemake@output[[2]] #"Stats/mg/mg.{read}_nonpareil.diversity.txt"

# If the input is empty just touch the output files
if (file.info(input_data)$size == 0) {
  print("Input file is empty. Touching output files.")
  system(paste("touch", outplot))
  system(paste("touch", outtab))
  q()
}

###################################################################################################
## Run curve 
###################################################################################################
print("Nonpareil calculations and curve.")
dir.create("Visualization", showWarnings = FALSE)
png(outplot,width=8.5,height=7,units = "cm",res=300,pointsize=7)
np_stats <- Nonpareil.curve(input_data)
dev.off()
print("Write Output")
write.table(data.frame("name"=c("model_a","model_b","nopareil_diversity"), 
                       "value"=c(summary(np_stats@model)$parameters[,1],np_stats$diversity),
                       stringsAsFactors=F),
             outtab,
             row.names=F,sep="\t",quote=F,col.names=F)

