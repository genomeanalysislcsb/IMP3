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

###################################################################################################
## Files from snakemake
###################################################################################################

cluster_dirs        <- snakemake@input #"Binning/selected_DASTool_bins/{i}/GRiD"
out_file            <- snakemake@output #"Stats/%s/bins.tsv" %ASS
binners             <- unlist(strsplit(snakemake@config[["binning"]][["binners"]],split=" "))

###################################################################################################
## Read DAStool output for all bins
###################################################################################################
print("Reading results for all binners.")
for(binner in binners){
  eval_file <- paste0("Binning/selected_",binner,".eval")
  if(file.exists(eval_file)){
    if(!"binTab" %in% ls()){
      binTab <- data.frame("binner"=binner,
                           read.delim(eval_file,stringsAsFactors=F),
                           stringsAsFactors = F)
    }else{
      binTab <- rbind(binTab,
                      data.frame("binner"=binner,
                                 read.delim(eval_file,stringsAsFactors=F),
                                 stringsAsFactors = F))
    }
  }
}

###################################################################################################
## Read DAStool output for selected bins
###################################################################################################
print("Reading DAStool results.")
selTab <- read.delim("Binning/selected_DASTool_summary.txt",stringsAsFactors=F)
binTab$selected_by_DASTool <- sapply(binTab$bin, function(x){
  if(x %in% selTab$bin){
    y <- x
  }else if(paste0(x,"_sub") %in% selTab$bin){
      y <- paste0(x,"_sub")
  }else{
      y <- ""
    }
  })

binTab <- merge(binTab,selTab,by.x="selected_by_DASTool",by.y="bin",all=T,suffix=c(".binner",".DASTool"))

###################################################################################################
## Read GRiD output for selected bins
###################################################################################################
print("Reading GRiD results.")
for(bin in setdiff(binTab$selected_by_DASTool[!is.na(binTab$selected_by_DASTool)],"")){
  if(!"gridTab" %in% ls()){
    gridTab <- data.frame("bin"=bin,
                         read.delim(paste0("Binning/selected_DASTool_bins/",bin,"/grid/reads.GRiD.txt"),
                                    stringsAsFactors=F)[,-1],
                         stringsAsFactors = F)
  }else{
    gridTab <- rbind(gridTab,
                    data.frame("bin"=bin,
                               read.delim(paste0("Binning/selected_DASTool_bins/",bin,"/grid/reads.GRiD.txt"),
                                          stringsAsFactors=F)[,-1],
                               stringsAsFactors = F))
  }
}

binTab <- merge(binTab,gridTab,by.x="selected_by_DASTool",by.y="bin",all=T)

###################################################################################################
## Save concatenated output
###################################################################################################

print("Writing concatenated bin table.")
write.table(binTab,out_file,sep="\t",quote=F,row.names=F)



