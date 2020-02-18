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
## Files from file system
###################################################################################################

ass                 <- snakemake@params[["ass"]] 
types               <- snakemake@params[["types"]] 
out_file_tsv        <- snakemake@output[[1]] #"Stats/all_stats.tsv"
out_file_rdat       <- snakemake@output[[1]] #"Stats/all_stats.Rdata"

readStat_files <- list()
for(type in types){
  readStats_files[[type]] <- paste0("Stats/",type,"/",type,".read_counts.txt")
}

###################################################################################################
## Read Preprocessing stats
###################################################################################################
print("Reading read stats.")
for(type in types){
  if(!"readTab" %in% ls()){
    readTab <- data.frame("type"=type,
                          read.delim(readStats_files[[type]],stringsAsFactors=F,header = F),
                          stringsAsFactors = F)
  }else{
    readTab <- rbind(readTab,
                     data.frame("type"=type,
                                read.delim(readStats_files[[type]],stringsAsFactors=F,header = F),
                                stringsAsFactors = F))
  }
}


###################################################################################################
## Save concatenated output
###################################################################################################

print("Writing concatenated read table.")
write.table(readTab,out_file_tsv,sep="\t",quote=F,row.names=F)
saveRDS(readTab,out_file_rdat,sep="\t",quote=F,row.names=F)


