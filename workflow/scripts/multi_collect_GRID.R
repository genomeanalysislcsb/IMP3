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

cluster_files       <- snakemake@input #"Dereplication/GRiD_multi/{i}/grid/reads.GRiD.txt"
out_file            <- snakemake@output #"Dereplication/allSamples_bin_stats.tsv"
binners             <- unlist(strsplit(snakemake@config[["dereplication"]][["cross_mapping_rebinning"]][["binners"]],split=" "))
samples             <- snakemake@params[["sample"]]

###################################################################################################
## Read DAStool output for all bins
###################################################################################################
print("Reading results for all binners.")
for(binner in binners){
  for(sample in samples){
    eval_file <- paste0("Dereplication/",sample,"/selected_",binner,".eval")
  }
  if(file.exists(eval_file)){
    if(!"binTab" %in% ls()){
      binTab <- data.frame("binner"=binner,
                           "sample"=sample,
                           read.delim(eval_file,stringsAsFactors=F),
                           stringsAsFactors = F)
    }else{
      binTab <- rbind(binTab,
                      data.frame("binner"=binner,
                                 "sample"=sample,
                                 read.delim(eval_file,stringsAsFactors=F),
                                 stringsAsFactors = F))
    }
  }
}
binTab$multiBin <- apply(binTab[,c("sample","bin")],1,function(x) paste(x,sep=".",collapse="."))

###################################################################################################
## Read DAStool output for selected bins
###################################################################################################
print("Reading DAStool results.")
for(sample in samples){
  eval_file <- paste0("Dereplication/",sample,"/selected_DASTool_summary.txt")
  if(file.exists(eval_file)){
    if(!"selTab" %in% ls()){
      selTab <- data.frame("sample"=sample,
                           read.delim(eval_file,stringsAsFactors=F))
    }else{
      selTab <- rbind(selTab,
                      data.frame("sample"=sample,
                                 read.delim(eval_file,stringsAsFactors=F)))
    }
  }}
selTab$multiBin <- apply(selTab[,c("sample","bin")],1,function(x) paste(x,sep=".",collapse="."))

binTab$selected_by_DASTool <- sapply(binTab$multiBin, function(x){
  if(x %in% selTab$multibin){
    y <- x
  }else if(paste0(x,"_sub") %in% selTab$multiBin){
      y <- paste0(x,"_sub")
  }else{
      y <- ""
    }
  })

binTab <- merge(binTab,selTab,by.x="selected_by_DASTool",by.y="multiBin",all=T,suffix=c(".binner",".DASTool"))

###################################################################################################
## Read GRiD output for selected bins
###################################################################################################
print("Reading GRiD results.")
for(bin in setdiff(binTab$selected_by_DASTool[!is.na(binTab$selected_by_DASTool)],"")){
  eval_file <- paste0("Dereplication/GRiD_multi/",bin,"/grid/reads.GRiD.txt")
  if(file.exists(eval_file)){
    if(!"gridTab" %in% ls()){
      gridTab <- data.frame("bin"=bin,
                            read.delim(eval_file,stringsAsFactors=F)[,-1],
                            stringsAsFactors = F)
    }else{
      gridTab <- rbind(gridTab,
                       data.frame("bin"=bin,
                                  read.delim(eval_file,stringsAsFactors=F)[,-1],
                                  stringsAsFactors = F))
    }
  }
}

binTab <- merge(binTab,gridTab,by.x="selected_by_DASTool",by.y="bin",all=T)

###################################################################################################
## Save concatenated output
###################################################################################################

print("Writing concatenated bin table.")
write.table(binTab,out_file[[1]],sep="\t",quote=F,row.names=F)



