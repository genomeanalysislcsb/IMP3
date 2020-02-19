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
## Functions
###################################################################################################

N50 <- function(lengths){
  tot <- sum(lengths)
  contls <- sort(lengths,decreasing = T)
  n50c <- max(which(cumsum(contls)<tot/2))
  n50 <- contls[n50c]
}

###################################################################################################
## Files from file system
###################################################################################################

steps               <- snakemake@params[["steps"]] 
ass                 <- snakemake@params[["ass"]] 
types               <- snakemake@params[["types"]] 
dbs                 <- unlist(strsplit(snakemake@params[["dbs"]],split=" "))
out_file_tsv        <- snakemake@output[[1]] #"Stats/all_stats.tsv"
out_file_rdat       <- snakemake@output[[2]] #"Stats/all_stats.Rdata"

###################################################################################################
## Read Preprocessing stats
###################################################################################################
if(file.exists(paste0("Stats/",type,"/",types[1],".read_counts.txt"))){
  readStats_files <- list()
  for(type in types){
    readStats_files[[type]] <- paste0("Stats/",type,"/",type,".read_counts.txt")
  }
  
  print("Reading read stats.")
  for(type in types){
    if(!"statsTab" %in% ls()){
      statsTab <- data.frame("type"=type,
                             read.delim(readStats_files[[type]],stringsAsFactors=F,header = F),
                             stringsAsFactors = F)
    }else{
      statsTab <- rbind(statsTab,
                        data.frame("type"=type,
                                   read.delim(readStats_files[[type]],stringsAsFactors=F,header = F),
                                   stringsAsFactors = F))
    }
  }
  colnames(statsTab)[2:3] <- c("step","number")
}
###################################################################################################
## Reads analysed by kraken
###################################################################################################
if(file.exists(paste0("Analysis/taxonomy/kraken/",types[1],".kraken.report"))){
  krakStats_files <- list()
  for(type in types){
    krakStats_files[[type]] <- paste0("Analysis/taxonomy/kraken/",type,".kraken.report")
  }
  
  print("Reading Kraken results.")
  for(type in types){
    krakRep <- read.delim(krakStats_files[[type]],stringsAsFactors=F,header = F)
    krakTemp <- data.frame("type"=type,
                           "step"=c("kraken classified reads","kraken unclassified reads",
                                    "kraken species","kraken genera"),
                           "number"=c(krakRep$V2[2],krakRep$V2[1],
                                      length(which(krakRep$V4=="S")),length(which(krakRep$V4=="G"))),
                           stringsAsFactors=F)
    rm(krakRep)
    if(!"statsTab" %in% ls()){
      statsTab <- krakTemp
    }else{
      statsTab <- rbind(statsTab,
                        krakTemp)
    }
  }
}
###################################################################################################
## Reads analysed by mOTUs
###################################################################################################
if(file.exists(paste0("Analysis/taxonomy/mOTUs/",types[1],".mOTU.counts.tsv"))){
  motuStats_files <- list()
  for(type in types){
    motuStats_files[[type]] <- paste0("Analysis/taxonomy/mOTUs/",type,".mOTU.counts.tsv")
  }
  
  print("Reading mOTU results.")
  for(type in types){
    motuRep <- read.delim(motuStats_files[[type]],stringsAsFactors=F,header = F, comment.char = "#")
    motuTemp <- data.frame("type"=type,
                           "step"=c("mOTU mapping reads",
                                    "mOTU species","mOTU genera"),
                           "number"=c(motuRep$V2[1],
                                      length(grep("s__",motuRep$V1)),length(grep("g__",motuRep$V1))),
                           stringsAsFactors=F)
    rm(motuRep)
    if(!"statsTab" %in% ls()){
      statsTab <- motuTemp
    }else{
      statsTab <- rbind(statsTab,
                        motuTemp)
    }
  }
}
###################################################################################################
## Reads mapping to Assembly
###################################################################################################
if(file.exists(paste0("Stats/",types[1],"/",ass,".assembly.contig_flagstat.txt"))){
  mapStats_files <- list()
  for(type in types){
    mapStats_files[[type]] <- paste0("Stats/",type,"/",ass,".assembly.contig_flagstat.txt")
  }
  
  #example:
  # 798852 + 0 in total (QC-passed reads + QC-failed reads) #lines in sam
  # 64233 + 0 secondary # 0x100 bit set = secondary alignment (all 0x4)
  # 0 + 0 supplementary
  # 0 + 0 duplicates
  # 621323 + 0 mapped (77.78% : N/A) # 0x4 bit not set = segment is mapping
  # 723816 + 0 paired in sequencing #input data
  # 361908 + 0 read1 #input data 
  # 361908 + 0 read2 #input data
  # 459758 + 0 properly paired (63.52% : N/A) # both 0x1 and 0x2 bits set and 0x4 bit not set = paired, each segment properly aligned according to the aligner
  # 539140 + 0 with itself and mate mapped #0x1 bit set and neither 0x4 nor 0x8 bits set = both segments of a pair are mapping somewhere 
  # 10955 + 0 singletons (1.51% : N/A) #both 0x1 and 0x8 bits set and bit 0x4 not set = paired, segment maps but partner doesn't
  # 76882 + 0 with mate mapped to a different chr 
  # 33150 + 0 with mate mapped to a different chr (mapQ>=5)
  
  print("Reading map stats.")
  for(type in types){
    mapTab <- read.delim(mapStats_files[[type]],stringsAsFactors=F,header = F)
    mapTemp <- data.frame("type"=type,
                          "step"=c("paired reads","se reads",
                                   "primary lines in bam",
                                   "segments mapped",
                                   "both partners map","one partner maps"),
                          "number"=c(mapTab$V1[7],mapTab$V1[1]-mapTab$V1[2]-mapTab$V1[6],
                                     mapTab$V1[1]-mapTab$V1[2],
                                     mapTab$V1[5],
                                     mapTab$V1[10]/2,mapTab$V1[11]),
                          stringsAsFactors=F)
    if(!"statsTab" %in% ls()){
      statsTab <- mapTemp
    }else{
      statsTab <- rbind(statsTab,
                        mapTemp)
    }
  }
}
###################################################################################################
## Reads mapping to Annotations
###################################################################################################
if(file.exists(paste0("Analysis/annotation/",types[1],".",dbs[1],".counts.tsv.summary"))){
  fcStats_files <- list()
  for(type in types){
    fcStats_files[[type]] <- list()
    for(db in dbs){
      fcStats_files[[type]][[db]] <- paste0("Analysis/annotation/",type,".",db,".counts.tsv.summary")
    }
  }
  
  print("Reading map stats for features.")
  for(type in types){
    mapTemp <- data.frame("type"=type,
                          "step"=paste("mapping to", dbs),
                          "number"=sapply(dbs,function(x){
                            mapTab <- read.delim(fcStats_files[[type]][[x]],stringsAsFactors=F,header = T)
                            mapTab[1,2]
                          }),
                          stringsAsFactors=F)
    if(!"statsTab" %in% ls()){
      statsTab <- mapTemp
    }else{
      statsTab <- rbind(statsTab,
                        mapTemp)
    }
  }
}
###################################################################################################
## Assembly stats
###################################################################################################
if(file.exists(paste0("Stats/",ass,"/",ass,".assembly.length.txt"))){
  print("Reading assembly stats.")
  contig_length <- read.delim(paste0("Stats/",ass,"/",ass,".assembly.length.txt"),
                              stringsAsFactors=F,header = F)
  colnames(contig_length) <- c("contig","length")
  contig_gc <- read.delim(paste0("Stats/",ass,"/",ass,".assembly.gc_content.txt"),
                          stringsAsFactors=F,header = T)
  colnames(contig_gc) <- c("contig","GC")
  
  depthAssStats_files <- list()
  breadthAssStats_files <- list()
  for(type in types){
    depthAssStats_files[[type]] <- paste0("Stats/",type,"/",ass,".assembly.contig_depth.txt")
    breadthAssStats_files[[type]] <- paste0("Stats/",type,"/",ass,".assembly.contig_coverage.txt")
  }
  
  assTemp <- data.frame("type"=rep(ass,4),
                        "step"=c("contigs","total assembly length",
                                 "N50",
                                 "mean contig length"),
                        "number"=c(nrow(contig_length),sum(contig_length$length),
                                   N50(contig_length$length),
                                   round(mean(contig_length$length))),
                        stringsAsFactors=F)
  
  if(!"statsTab" %in% ls()){
    statsTab <- assTemp
  }else{
    statsTab <- rbind(statsTab,
                      assTemp)
  }
  
  depthAssStat <- list()
  breadthAssStat <- list()
  for(type in types){
    depthAssStat[[type]] <- read.delim(depthAssStats_files[[type]],stringsAsFactors=F,header = F)
    breadthAssStat[[type]] <- read.delim(breadthAssStats_files[[type]],stringsAsFactors=F,header = F)[,c(1,7)]
    colnames(depthAssStat[[type]]) <- c("contig","depth")
    colnames(breadthAssStat[[type]]) <- c("contig","breadth")
    covTemp <- data.frame("type"=type,
                          "step"=c("contigs_max_depth",
                                   "contigs_mean_depth>0","contigs_mean_depth>1",
                                   "contigs_breadth>0","contigs_breadth=1"),
                          "number"=c(round(max(depthAssStat[[type]]$depth)),
                                     nrow(depthAssStat[[type]]),length(which(depthAssStat[[type]]$depth>=1)),
                                     length(which(breadthAssStat[[type]]$breadth>0)),length(which(breadthAssStat[[type]]$breadth==1))),
                          stringsAsFactors=F)
    if(!"statsTab" %in% ls()){
      statsTab <- covTemp
    }else{
      statsTab <- rbind(statsTab,
                        covTemp)
    }
  }
}
###################################################################################################
## Contigs analysed by kraken
###################################################################################################
if(file.exists(paste0("Analysis/taxonomy/kraken/",ass,".kraken.parsed.tsv"))){
  print("Reading Kraken contig results.")
  krakAssFile <- paste0("Analysis/taxonomy/kraken/",ass,".kraken.parsed.tsv")
  krakTemp <- data.frame("type"=ass,
                         "step"=c("kraken classified contigs"),
                         "number"=c(as.numeric(unlist(strsplit(system2("wc",
                                                                       args=c("-l",krakAssFile),
                                                                       stdout=T),split=" "))[1])-1),
                         stringsAsFactors=F)
  if(!"statsTab" %in% ls()){
    statsTab <- krakTemp
  }else{
    statsTab <- rbind(statsTab,
                      krakTemp)
  }
}
###################################################################################################
## Annotation stats
###################################################################################################
if(file.exists(paste0("Stats/",types[1],"/annotation/",types[1],".gene.len"))){
  print("Reading annotation stats.")
  gene_length <- read.delim(paste0("Stats/",types[1],"/annotation/",types[1],".gene.len"),
                            stringsAsFactors=F,header = F)
  colnames(gene_length) <- c("gene","length")
  
  if(file.exists("Analysis/annotation/annotation_CDS_RNA_hmms.gff")){
    gene_gff <- read.delim("Analysis/annotation/annotation_CDS_RNA_hmms.gff",
                           stringsAsFactors=F,header = F)
  }else{
    gene_gff <- read.delim("Analysis/annotation/annotation.filt.gff",
                           stringsAsFactors=F,header = F)
  }
  colnames(gene_gff) <- c("contig","source","feature","start","end","score","strand","frame","attribute")
  
  depthGeneStats_files <- list()
  for(type in types){
    depthGeneStats_files[[type]] <- paste0("Stats/",type,"/annotation/",type,".gene_depth.avg")
  }
  
  featureCount_files <- list()
  for(db in dbs){
    featureCount_files[[db]] <- paste0("Analysis/annotation/",type[1],".",db,".counts.tsv")
  }
  
  annTemp <- data.frame("type"=ass,
                        "step"=c("genes","mean gene length",
                                 "CDSs","rRNAs","tmRNAs","tRNAs",
                                 dbs, paste(dbs,"categories")),
                        "number"=c(nrow(gene_length),round(mean(gene_length$length)),
                                   length(which(gene_gff$feature=="CDS")),
                                   length(which(gene_gff$feature=="rRNA")),
                                   length(which(gene_gff$feature=="tmRNA")),
                                   length(which(gene_gff$feature=="tRNA")),
                                   sapply(dbs,function(x){
                                     patt=paste0(";",x,"=")
                                     length(grep(patt,gene_gff))
                                   }),
                                   sapply(dbs,function(x){
                                     fileoi <- featureCount_files[[x]]
                                     as.numeric(unlist(strsplit(system2("wc",
                                                                        args=c("-l",fileoi),
                                                                        stdout=T),split=" "))[1])-1
                                   })),
                        stringsAsFactors=F)
  
  if(!"statsTab" %in% ls()){
    statsTab <- annTemp
  }else{
    statsTab <- rbind(statsTab,
                      annTemp)
  }
  
  depthGeneStat <- list()
  for(type in types){
    depthGeneStat[[type]] <- read.delim(depthGeneStats_files[[type]],stringsAsFactors=F,header = F)
    colnames(depthGeneStat[[type]]) <- c("gene","depth")
    covTemp <- data.frame("type"=type,
                          "step"=c("max_gene_depth","genes_mean_depth>0"),
                          "number"=c(round(max(depthGeneStat[[type]]$depth)),length(which(depthGeneStat[[type]]$depth>0))),
                          stringsAsFactors=F)
    if(!"statsTab" %in% ls()){
      statsTab <- covTemp
    }else{
      statsTab <- rbind(statsTab,
                        covTemp)
    }
  }
}
###################################################################################################
## Binning stats
###################################################################################################
if(file.exists(paste0("Stats/",ass,"/",ass,".bins.tsv"))){
  print("Reading binning stats.")
  binTab <- read.delim(paste0("Stats/",ass,"/",ass,".bins.tsv"),
                       stringsAsFactors=F,header = T)
  
  rnas_bin_tab <- matrix(0,nrow=length(binTab$selected_by_DASTool[binTab$selected_by_DASTool!=""]),
                         ncol=6,
                         dimnames=list(binTab$selected_by_DASTool[binTab$selected_by_DASTool!=""],
                                       c("complete_16S","partial_16S","complete_23S","partial_23S",
                                         "any_5S","tRNAs")))
  for(bin in binTab$selected_by_DASTool[binTab$selected_by_DASTool!=""]){
    contigs_temp <- read.delim(paste0("Binning/selected_DASTool_bins/",bin,"/contigs.bed"),
                               stringsAsFactors=F,header = F)[,1]
    products_bin <- gsub(".+;product=","",gene_gff[gene_gff$contig %in% contigs_temp & gene_gff$feature != "CDS",9])
    if(length(products_bin)>0){
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,1] <- length(which(grepl("16S ribosomal RNA",products_bin) &
                                                                    !grepl("partial",products_bin)))
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,2] <- length(which(grepl("16S ribosomal RNA",products_bin) &
                                                                    grepl("partial",products_bin)))
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,3] <- length(which(grepl("23S ribosomal RNA",products_bin) &
                                                                    !grepl("partial",products_bin)))
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,4] <- length(which(grepl("23S ribosomal RNA",products_bin) &
                                                                    grepl("partial",products_bin)))
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,5] <- length(which(grepl("5S ribosomal RNA",products_bin)))
      rnas_bin_tab[rownames(rnas_bin_tab)==bin,6] <- length(unique(grep("tRNA",products_bin,value = T)))
    }
    rm(contigs_temp)
    rm(products_bin)
  }
  
  binTab <- merge(binTab,rnas_bin_tab,by.x="selected_by_DASTool",by.y=0,all.x=T)
  
  binTemp <- data.frame("type"=ass,
                        "step"=c("total_bins","DASTool_selected_bins",
                                 "HQ_drafts",
                                 "HQ_drafts_partial_RNAs",
                                 "HQ_drafts_minus_RNAs",
                                 "MQ_drafts",
                                 "LQ_drafts",
                                 "binned_length_DASTool",
                                 "binned_contigs_DASTool"),
                        "number"=c(nrow(binTab),length(which(binTab$selected_by_DASTool!="")),
                                   length(which(binTab$SCG_completeness.DASTool>90 &
                                                  binTab$SCG_redundancy.DASTool<5 &
                                                  binTab$complete_16S>0 &
                                                  binTab$complete_23S>0 &
                                                  binTab$any_5S>0 &
                                                  binTab$tRNAs >= 18)),
                                   length(which(binTab$SCG_completeness.DASTool>90 &
                                                  binTab$SCG_redundancy.DASTool<5 &
                                                  (binTab$complete_16S+binTab$partial_16S)>0 &
                                                  (binTab$complete_23S+binTab$partial_23S)>0 &
                                                  !(binTab$complete_16S>0 &
                                                      binTab$complete_23S>0) &
                                                  binTab$any_5S>0 &
                                                  binTab$tRNAs >= 18)),
                                   length(which(binTab$SCG_completeness.DASTool>90 &
                                                  binTab$SCG_redundancy.DASTool<5 & 
                                                  !((binTab$complete_16S+binTab$partial_16S)>0 &
                                                      (binTab$complete_23S+binTab$partial_23S)>0 &
                                                      binTab$any_5S>0 &
                                                      binTab$tRNAs >= 18))),
                                   length(which(binTab$SCG_completeness.DASTool>=50 &
                                                  binTab$SCG_redundancy.DASTool<10 & 
                                                  !(binTab$SCG_completeness.DASTool>90 &
                                                      binTab$SCG_redundancy.DASTool<5))),
                                   length(which(binTab$SCG_completeness.DASTool<50 &
                                                  binTab$SCG_redundancy.DASTool<10 )),
                                   sum(binTab$size.DASTool[binTab$selected_by_DASTool!=""]),
                                   sum(binTab$contigs.DASTool[binTab$selected_by_DASTool!=""])),
                        stringsAsFactors=F)
  
  if(!"statsTab" %in% ls()){
    statsTab <- binTemp
  }else{
    statsTab <- rbind(statsTab,
                      binTemp)
  }
}
###################################################################################################
## Save concatenated output
###################################################################################################

print("Writing concatenated read table.")
write.table(statsTab,out_file_tsv,sep="\t",quote=F,row.names=F)
save.image(out_file_rdat)

