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

cluster_files       <- snakemake@input[["dRep"]] #"Chdb.csv","Cdb.csv","Sdb.csv","Wdb.csv","Widb.csv","Ndb.csv"
samples             <- snakemake@params[["sample"]] #sample names
stats_files         <- snakemake@input[["ori"]] # "all_bin_stats.tsv" of original samples
out_files            <- snakemake@output #"dereplicated_bin_stats.tsv","originalNdereplicated_bin_stats.tsv"

###################################################################################################
## Read dRep output for all bins
###################################################################################################

bins <- read.delim(cluster_files[[1]],sep=",",stringsAsFactors = F)
bins <- merge(bins,read.delim(cluster_files[[2]],sep=",",stringsAsFactors = F),by=1,all=T)
bins <- merge(bins,read.delim(cluster_files[[3]],sep=",",stringsAsFactors = F),by=1,all=T)
wins <- read.delim(cluster_files[[4]],sep=",",stringsAsFactors = F)
bins$dRepWinner <- sapply(bins$Bin.Id,function(x) if(x %in% wins$genome) x else NA)
bins <- merge(read.delim(cluster_files[[5]],sep=",",
                   stringsAsFactors = F)[,c("genome","cluster_members",
                                            "closest_cluster_member","furthest_cluster_member")],
              bins, by=1,all=T)
wins <- bins[!is.na(bins$dRepWinner),]
wins <- wins[,-ncol(wins)]
wins$memberNames <- sapply(wins$secondary_cluster,
                           function(x) {
                             mem <- setdiff(bins$genome[bins$secondary_cluster==x],
                                            wins$genome)
                             mem <- mem[!is.na(mem)]
                             if(length(mem)>0) r <- paste(mem,sep=";",collapse=";") else r<- ""
                             r})
dist <- read.delim(cluster_files[[5]],sep=",",stringsAsFactors = F)
bins$distToWinner <- sapply(bins$genome,function(x){
  if(x %in% wins$genome) r <- 0 else{
    if(! x %in% dist$querry) r <- NA else{
      r <- dist$ani[dist$querry == x & 
                      dist$reference == wins$genome[wins$secondary_cluster==bins$secondary_cluster[bins$genome==x]]]
    }
  }
r
})
wins$genome <- gsub(".contigs.fa","",wins$genome)
bins$genome <- gsub(".contigs.fa","",bins$genome)
bins$dRepWinner <- gsub(".contigs.fa","",bins$dRepWinner)

###################################################################################################
## Merge with previous data
###################################################################################################
for(s in 1:length(samples)){
  stat <- read.delim(stats_files[[s]],stringsAsFactors = F)
  stat <- stat[!is.na(stat$selected_by_DASTool) & stat$selected_by_DASTool != "",]
  if(nrow(stat)>0){
    stat$selected_by_DASTool <- paste(samples[[s]],stat$selected_by_DASTool,sep=".")
    if("allS" %in% ls()){
      allS <- rbind(allS,stat)
    }else{
      allS <- stat
    }
  }
}
bins <- merge(bins,allS,by=1,all.x=T,suffixes=c(".dRep",".ori"))
wins <- merge(wins,allS,by=1,all.x=T,suffixes=c(".dRep",".ori"))

###################################################################################################
## Save output
###################################################################################################

write.table(bins,out_files[[1]],sep="\t",quote=F,row.names=F)
write.table(wins,out_files[[2]],sep="\t",quote=F,row.names=F)

