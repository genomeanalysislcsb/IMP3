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
library(alluvial)
library(data.table)

###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input[["stats"]] #"Stats/all_stats.Rdata"
bin_memberships     <- snakemake@input[["bins"]]  # Binning/*/scaffold2bin.tsv
mag_memberships     <- snakemake@input[["dastool"]]  # Binning/selected_DASTool_scaffolds2bin.tsv
outplot             <- snakemake@output # "Visualization/Assembly_length_hist.png"
#                                         Visualization/Assembly_cov_hist.png"
                                        #"Visualization/Assembly_gc_hist.png"
                                         #"Visualization/Assembly_cov_gc_scatter.png"

###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

binners <- unique(binTab$binner)

###################################################################################################
## plot
###################################################################################################
print("Plotting alluvial diagram.")
dir.create("Visualization", showWarnings = FALSE)

plotfile <- outplot[[1]]
cont2bin <- list()
for(binner in binners){
  cont2bin[[binner]] <- read.delim(grep(binner,bin_memberships,value=T),stringsAsFactors = F,header=F)
  colnames(cont2bin[[binner]]) <- c("contig","bin")
}
cont2bin[["DASTool"]] <- read.delim(mag_memberships,stringsAsFactors = F,header=F)
colnames(cont2bin[["DASTool"]]) <- c("contig","bin")

cont2bin <- rbindlist(lapply(1:length(cont2bin),function(x) data.frame(cont2bin[[x]],
                                                                       "tool"=names(cont2bin)[x],
                                                                       stringsAsFactors = F)))
tools <- levels(as.factor(cont2bin$tool))
cont2bin <- apply(tapply(cont2bin$bin,list(cont2bin$contig,cont2bin$tool),paste),1,function(x)paste(x,sep="|",collapse="|"))
cont2bin <- merge(cont2bin,contig_length,by.x=0,by.y=1,all.x=T)
cont2bin <- aggregate(cont2bin$length,list(cont2bin$x),sum)
cont2bin$Group.1 <- as.character(cont2bin$Group.1)
cont2bin <- data.frame(t(sapply(cont2bin$Group.1,function(x) unlist(strsplit(x,split="|",fixed=T)))),
                       "freq"=cont2bin$x,stringsAsFactors = F)
colnames(cont2bin)[1:length(tools)] <- tools
if("MaxBin" %in% colnames(cont2bin)){
  cont2bin$MaxBin <- gsub("maxbin_res.","",gsub(".fasta","",cont2bin$MaxBin))
  cont2bin$DASTool <- gsub("maxbin_res.","",gsub(".fasta","",cont2bin$DASTool))
}

png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
par("mar"=c(0.2,0.2,0.2,0.2),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
alluvial(cont2bin[,c(binners,"DASTool")],freq = cont2bin$freq,
         col=colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(cont2bin$DASTool)))[as.numeric(as.factor(cont2bin$DASTool))],
         border=NA)
dev.off()
