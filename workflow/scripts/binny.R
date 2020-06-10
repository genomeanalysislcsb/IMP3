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
print("Loading required R libraries")
library(genomeIntervals)
library(ggplot2)
library(gtools)
library(gtable)
library(data.table)
library(reshape)
library(grid)
library(grDevices)
library(stringr)
library(xtable)
#binny libraries:
library(caTools)
library(fpc)
library(FNN)
library(RColorBrewer)
library(scales)
library(diptest)
library(mixtools)
library(gclus)

###################################################################################################
## Files now come from snakemake
###################################################################################################

cluster_dir         <- snakemake@input[["outdir"]]
binny_dir           <- snakemake@params[["binnydir"]]
MG.depth_file       <- snakemake@input[["mgdepth"]]
coords_file         <- snakemake@input[["vizbin"]]
annot_file          <- snakemake@input[["gff"]]
function_script     <- snakemake@params[["plot_functions"]]
pk 		    <- snakemake@config[["binning"]][["binny"]][["pk"]]
nn 		    <- snakemake@config[["binning"]][["binny"]][["nn"]]
all_outputs <- snakemake@output

###################################################################################################
## Initialize functions for various calculations and normalizations
###################################################################################################
print("START: Reading functions")
source(function_script)
print("DONE: Reading functions")

###################################################################################################
## Read in the necessary input files
###################################################################################################
print("START: Reading data")

## MG depth (bedtools coverageBed output)
print("Read in MG depth file")
MG.depth <- read.table(MG.depth_file, colClasses=c("factor", "numeric"))
colnames(MG.depth) <- c("contig", "MG_depth")

## gff annotation file
print("Read in gff3 annotation file")
annot <- readGff3(annot_file, isRightOpen = FALSE, quiet=TRUE)

# extract only interesting information
print("Processing gff3 annotation file")
dbs <- setdiff(unique(unlist(lapply(parseGffAttributes(annot),names))),c("ID","inference","locus_tag","product","note","Name","gene"))
annot.1 <- as.data.frame(
  cbind(
    as.character(annot@annotation$seq_name),
    annot@.Data,
    getGffAttribute(annot, dbs)
  )
)
colnames(annot.1)[1:3] <- c("contig", "start", "end")

# Remove rows with no start and stop positions
if (nrow(annot.1[is.na(annot.1$start),]) != 0) {
  annot.1 <- annot.1[-is.na(annot.1$start),]
} else if (nrow(annot.1[is.na(annot.1$end),]) != 0) {
  annot.1 <- annot.1[-is.na(annot.1$end),]
} else {
  print("Genes are OK!")
}

# Replace non-annotated genes with NA
annot.1[,2:3] <- sapply(annot.1[2:3], as.character)
annot.1[,2:3] <- sapply(annot.1[2:3], as.numeric)
save.image(paste0(binny_dir,"/binny_WS.Rdata"))

# create annotation table
if(!"essential" %in% colnames(annot.1)){
  for(o in setdiff(unlist(all_outputs),paste0(binny_dir,"/binny_WS.Rdata"))){
    system(paste0("touch ",o))
  }
  print("no essential genes. Ending Binny")
  system("touch binny.error")
}else{
  print("Creating annotation table")
  annot.E.prelim <- annot.1[,c("contig","start","end","essential")]
  annot.E.prelim <- annot.E.prelim[!is.na(annot.E.prelim$essential),]
  annot.E.prelim$essential <- as.character(annot.E.prelim$essential)
  annot.E.prelim$contig <- as.character(annot.E.prelim$contig)
  save.image(paste0(binny_dir,"/binny_WS.Rdata"))
  
  print("Reading in vizbin coordinates")
  coords <- read.table(coords_file, colClasses=c("factor", "numeric", "numeric"),
                       sep="\t", col.names=c("contig", "x", "y"))
  
  if(nrow(coords)==0 | length(intersect(annot.E.prelim$contig,coords$contig))==0 ){
    for(o in setdiff(unlist(all_outputs),paste0(binny_dir,"/binny_WS.Rdata"))){
      system(paste0("touch ",o))
    }
    system("touch binny.error")
    print("no vizbin coordinates for essential genes. Ending Binny")
  }else{
    print("DONE: Reading data")
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###################################################################################################
    ## Incorporate coordinates from vizbin
    ###################################################################################################
    print("START: Incorporating VizBin and annotation data")
    
    vb_dat <- merge(MG.depth, coords, by=c("contig"), all=T, incomparables=NA)
    vb_dat <- vb_dat[!is.na(vb_dat$x),]
    vb_dat$MG_depth[is.na(vb_dat$MG_depth)] <- 0
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    contigInfo <- vb_dat
    contigInfo$contig <- as.character(contigInfo$contig)
    annot.E <- annot.E.prelim[annot.E.prelim$contig %in% vb_dat$contig,]
    essAnno <- aggregate(annot.E$essential,list(annot.E$contig),function(x) paste(x,sep=";",collapse=";"))
    colnames(essAnno) <- c("contig","essentialGene")
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    contigInfo <- merge(contigInfo,essAnno,by.x=1,by.y=1,all.x=T)
    contigInfo$essentialGene[is.na(contigInfo$essentialGene)] <- "notEssential"
    
    print("DONE: Incorporating VizBin and annotation data")
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    
    ###################################################################################################
    ## binny
    ###################################################################################################
    
    print("START: Running binny on VizBin data")
    # Handle missing and infinite values
    #vb_dat[is.na(vb_dat)] <- 0
    
    # auto clustering
    coco <- c(which(colnames(contigInfo)=="x"),which(colnames(contigInfo)=="y"))
    #estimate reachability distance (based on number of neighbouting points)
    skn <- sort(knn.dist(contigInfo[,coco],pk)[,pk])
    sdkn <- runsd(skn,10)
    est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))]
    #document estimated reachability distances
    write.table(t(c("scan","reachabilityDist")),
                paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                row.names=F,col.names=F,quote=F,sep="\t")
    write.table(t(c("first",est)),
                paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                row.names=F,col.names=F,quote=F,sep="\t",append=T)
    #DBscan
    cdb <- dbscan(contigInfo[,coco],est,pk)
    
    #count contigs in clusters and document clustering steps
    pdf(paste0(binny_dir,paste("/scatterPlot1",pk,nn,"pdf",sep=".")))
    plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
    j <- 1
    cdbTab <- data.frame("cluster"=names(table(cdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
    for(i in names(table(cdb$cluster))) {
      points(contigInfo[cdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
      cdbTab$contigs[cdbTab$cluster==i] <- table(cdb$cluster)[names(table(cdb$cluster))==i]
      cdbTab$numberEss[cdbTab$cluster==i] <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
      cdbTab$uniqueEss[cdbTab$cluster==i] <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
      j<-j+1
    }
    box()
    dev.off()
    
    write.table(cdbTab,paste0(binny_dir,paste("/clusterFirstScan",pk,nn,"tsv", sep=".")),
                sep="\t",row.names=F,quote=F)
    write.table(t(c("clusterName","cutoff")),paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),
                sep="\t",row.names=F,col.names=F,quote=F)
    
    #assign clusters to contigs 
    ### names are "N" for contigs recognized as noise
    #### "E" for clusters without essential genes
    #### "C" for clusters with less than or equal to 20% duplicated essential genes
    #### "D" for clusters with more than 20% duplicated essential genes
    clusterRes <- data.frame("contig"=contigInfo$contig,"cluster"="x",stringsAsFactors=F)
    clusterRes$cluster[cdb$cluster==0] <- "N"
    
    emptyClus <- cdbTab$cluster[cdbTab$uniqueEss==0&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% emptyClus] <- paste("E",cdb$cluster[cdb$cluster %in% emptyClus],sep="")
    
    closeClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss<=1.2&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% closeClus] <- paste("C",cdb$cluster[cdb$cluster %in% closeClus],sep="")
    
    duClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss>1.2&cdbTab$cluster!=0]
    clusterRes$cluster[cdb$cluster %in% duClus] <- paste("D",cdb$cluster[cdb$cluster %in% duClus],sep="")
    
    #cut clusters by metagenomic coverage depth, document cut-offs
    write.table(t(c("cluster","cutoff")),
                paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),
                sep="\t",col.names=F,row.names=F,quote=F)
    for(ds in unique(clusterRes$cluster[grep("D",clusterRes$cluster)])){
      clusterRes <- muClus(ds,clusterRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1)
    }
    
    #second iteration of clustering on clusters with more than 20% duplicated essential genes 
    ### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
    ### - neighbouring points are increased by 2
    pk2 <- pk + 2
    pdf(paste0(binny_dir,paste("/scatterPlots2",pk,nn,"pdf",sep=".")))
    for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
      bbInfo <- contigInfo[clusterRes$cluster==bb,]
      if(nrow(bbInfo)>nn){
        sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
        if(length(sknBB>10)){
          sdknBB <- runsd(sknBB,10)
          plot(sdknBB,sknBB)
          plot(sdknBB)
          abline(h=quantile(sdknBB,0.975))
          estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
          if(is.na(estBB)){
            estBB <- sort(sknBB)[max(sdknBB)]
            print("estBB could not be found")
          }
          plot(sknBB)
          abline(h=estBB)
          
          write.table(t(c(bb,estBB)),
                      paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                      row.names=F,col.names=F,quote=F,sep="\t",append=T)
          BBcdb <- dbscan(bbInfo[,coco],estBB,pk2)
          
          plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
          j <- 1
          BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
          for(i in names(table(BBcdb$cluster))) {
            points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
            BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
            BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
            BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
            j<-j+1
          }
          write.table(BBcdbTab,
                      paste0(cluster_dir,paste("/blob1Cluster",pk,nn,gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                                                                  bb)),"tsv",sep=".")),
                      sep="\t",row.names=F,quote=F)
          
          BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
          BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
          clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
          
          BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
          BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
          clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                          gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          
          BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
          BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
          if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
          clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                          gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          
          BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
          BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
          uniD <- unique(paste(gsub("B","D",bb),
                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                               sep="."))
          clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                       gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                       sep=".")
          for(bds in uniD){
            clusterRes <- muClus(bds,clusterRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1)
          }
        }else{
          print(paste("Cluster",bb,"too small for knn-approach.\n"))
        }
      }else{
        print(paste("Cluster",bb,"too small for knn-approach.\n"))
      }
    }
    dev.off()
    
    #third iteration of clustering on clusters with more than 20% duplicated essential genes
    ### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
    ### - neighbouring points are increased by 2 again
    pk4 <- pk + 4
    pdf(paste0(binny_dir,paste("/scatterPlots3",pk,nn,"pdf",sep=".")))
    for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
      bbInfo <- contigInfo[clusterRes$cluster==bb,]
      if(nrow(bbInfo)>nn){
        sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
        if(length(sknBB>10)){
          sdknBB <- runsd(sknBB,10)
          plot(sdknBB,sknBB)
          plot(sdknBB)
          abline(h=quantile(sdknBB,0.975))
          estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
          if(is.na(estBB)){
            estBB <- sort(sknBB)[max(sdknBB)]
            print("estBB could not be found")
          }
          plot(sknBB)
          abline(h=estBB)
          
          write.table(t(c(bb,estBB)),
                      paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                      row.names=F,col.names=F,quote=F,sep="\t",append=T)
          BBcdb <- dbscan(bbInfo[,coco],estBB,pk4)
          
          plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
          j <- 1
          BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
          for(i in names(table(BBcdb$cluster))) {
            points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
            BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
            BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
            BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
            j<-j+1
          }
          write.table(BBcdbTab,
                      paste0(cluster_dir,paste("/blob2Cluster",pk,nn,gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                                                                  bb)),"tsv",sep=".")),
                      sep="\t",row.names=F,quote=F)
          
          BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
          BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
          clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
          
          BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
          BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
          clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                          gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          
          BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
          BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
          if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
          clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                          gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
          BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
          uniD <- unique(paste(gsub("B","D",bb),
                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                               sep="."))
          clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                       gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                       sep=".")
          for(bds in uniD){
            clusterRes <- muClus(bds,clusterRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1)
          }
        }else{
          print(paste("Cluster",bb,"too small for knn-approach.\n"))
        }
      }else{
        print(paste("Cluster",bb,"too small for knn-approach.\n"))
      }
    }
    dev.off()
    
    #fourth iteration of clustering on clusters with more than 20% duplicated essential genes
    ### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
    ### - neighbouring points are increased by 2 again
    pk6 <- pk + 6
    pdf(paste0(binny_dir,paste("/scatterPlots4",pk,nn,"pdf",sep=".")))
    for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
      bbInfo <- contigInfo[clusterRes$cluster==bb,]
      if(nrow(bbInfo)>nn){
        sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
        if(length(sknBB>10)){
          sdknBB <- runsd(sknBB,10)
          plot(sdknBB,sknBB)
          plot(sdknBB)
          abline(h=quantile(sdknBB,0.975))
          estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
          if(is.na(estBB)){
            estBB <- sort(sknBB)[max(sdknBB)]
            print("estBB could not be found")
          }
          plot(sknBB)
          abline(h=estBB)
          
          write.table(t(c(bb,estBB)),
                      paste0(binny_dir,paste("/reachabilityDistanceEstimates",pk,nn,"tsv",sep=".")),
                      row.names=F,col.names=F,quote=F,sep="\t",append=T)
          BBcdb <- dbscan(bbInfo[,coco],estBB,pk6)
          
          plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
          j <- 1
          BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
          for(i in names(table(BBcdb$cluster))) {
            points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
            BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
            BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
            BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
            j<-j+1
          }
          write.table(BBcdbTab,
                      paste0(cluster_dir,paste("/blob3Cluster",pk,nn,gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                                                                  bb)),"tsv",sep=".")),
                      sep="\t",row.names=F,quote=F)
          
          BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
          BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
          clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
          
          BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
          BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
          clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                          gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          
          BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
          BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
          if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
          clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                          gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                          sep=".")
          
          BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
          BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
          BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
          uniD <- unique(paste(gsub("B","D",bb),
                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                               sep="."))
          clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                       gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                       sep=".")
          for(bds in uniD){
            clusterRes <- muClus(bds,clusterRes,paste0(binny_dir,paste("/bimodalClusterCutoffs",pk,nn,"tsv",sep=".")),1)
          }
        }else{
          print(paste("Cluster",bb,"too small for knn-approach.\n"))
        }
      }else{
        print(paste("Cluster",bb,"too small for knn-approach.\n"))
      }
    }
    dev.off()
    
    #assignment of final names to clusters
    ### based on completeness
    #### "P": more than 100/109 essential genes, less than 115 essential genes in total (<14% duplicated genes, >92% complete)
    #### "G": more than 71/109 essential genes (>65% complete), less than 20% in duplicate
    #### "O": more than 51/109 essential genes (>47% complete), less than 20% in duplicate
    #### "L": more than 31/109 essential genes (>28% complete), less than 20% in duplicate
    #### "C": at least 1/109 essential genes (>= 1% complete), less than 20% in duplicate
    #### "E": no essential genes
    #### "B": at least 1/109 essential genes (>= 1% complete), at least 20% in duplicate
    #### "N": noise
    for(clus in grep("C",unique(clusterRes$cluster),value=T)){
      uni <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";")))))
      ess <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";"))))))
      if(uni>100&ess<115){
        clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","P",clus)
      } else if(uni>71){
        clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","G",clus)
      } else if(uni>51){
        clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","O",clus)
      } else if(uni>31){
        clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","L",clus)
      }
    }
    
    #save results
    clusterRes$cluster <- gsub("(2.){8,}","2222.",gsub("(1.){8,}","1111.",
                                                       clusterRes$cluster))
    write.table(clusterRes,
                paste0(binny_dir,paste("/contigs2clusters",pk,nn,"tsv",sep=".")),
                sep="\t",row.names=F,quote=F)
    # write.table(clusterRes[grep("^[PGOL]",clusterRes$cluster),],
    #            "scaffold2bin.tsv",
    #            sep="\t",row.names=F,quote=F)
    save(list=c("contigInfo","coco","pk","nn","essAnno","find.cutoff","muClus","clusterRes"),
         file=paste0(binny_dir,paste("/clusteringWS",pk,nn,"Rdata",sep=".")))
    saveRDS(clusterRes,paste0(binny_dir,paste("/contigs2clusters",pk,nn,"RDS",sep=".")))
    
    print("Saved tables. Plotting final map now.")
    #plot final map
    essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(111)[111:1]
    png(paste0(binny_dir,paste("/finalClusterMap",pk,nn,"png",sep=".")),width=17,height=17,pointsize=7,res=300,units = "cm")
    layout(matrix(1:2,nrow=1,ncol=2),c(5/6,1/6))
    par(mar=rep(1,4),pty="s")
    plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F,bty="o")
    cI <- data.frame("essNum"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                     "medX"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                     stringsAsFactors=F)
    i <- 1
    for(clus in unique(clusterRes$cluster)){
      uni <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";"))))))
      cI$essNum[i] <- uni
      rownames(cI)[i] <- clus
      cI$medX[i] <- median(contigInfo[clusterRes$cluster==clus,coco[1]])
      cI$medY[i] <- median(contigInfo[clusterRes$cluster==clus,coco[2]])
      cI$high[i] <- max(contigInfo[clusterRes$cluster==clus,coco[2]])
      rm(uni)
      i <- i +1
    }
    cI2 <- cI[rownames(cI)!="N",][order(cI$essNum[rownames(cI)!="N"]),]
    cI2 <- rbind(cI[rownames(cI)=="N",],cI2)
    for(i in 1:nrow(cI2)){
      clus <- rownames(cI2)[i]
      if(clus=="N") pc <- "grey" else if(grepl("E",clus)) pc <- "black" else if(grepl("B",clus)) pc <- brewer.pal(8,"Set3")[8] else{
        pc <- essPal[cI2$essNum[i]]
      }
      points(contigInfo[clusterRes$cluster==clus,coco],cex=0.3,pch=16,col=pc)
      if(cI2$essNum[i]>15&clus!="N") text(cI2$medX[i],cI2$high[i]+0.5,labels=clus,cex=1.0,font=2,col=pc,pos=3)
      rm(pc)
    }
    par(mar=c(3,3,3,1),pty="m")
    plotcolors(cbind(c(colorRampPalette(brewer.pal(11,"Spectral"))(109),"black"),c(colorRampPalette(brewer.pal(11,"Spectral"))(109),"black")))
    hi <- hist(cI$essNum,breaks=0:109,plot=F)
    a <- apply(cbind(hi$counts,hi$breaks[-length(hi$breaks)]+1),1,function(x) lines(c(0.5,(x[1]*1.9/max(hi$counts[-1]))+0.5),c(rep(x[2],2))))
    abline(h=c(101.5,71.5,51.5,31.5,1.5),col=c(rep("grey20",4),"grey80"),lty=2)
    axis(1,seq(from=0.5,to=2.4,length.out=3),labels=c(0,max(hi$counts[-1])/2,max(hi$counts[-1])))
    axis(2,seq(from=1,to=11000/100,length.out=6),labels=c(0,20,40,60,80,100),las=1)
    mtext("% completeness",2,2)
    mtext("clusters",1,2)
    dev.off()
    print("binny done")
    
    rm(list=c("contigInfo","coco","pk","nn","clusterRes","essPal","cI","cI2","hi","a","clus","ess","BBduClus","BBclusterRes","BBduCont","uniD","BBcloseClus","BBcloseCont","BBemptyClus","BBemptyCont","BBcdbTab","bbInfo","estBB","BBcdb","sknBB","emptyClus","closeClus","cdbTab","skn","sdkn","est","sdknBB","pk2","pk4","pk6"))
    
    
    ####################################################################
    ## Save the R workspace
    ####################################################################
    
    print("START: Saving R image: binny_WS.Rdata")
    save.image(paste0(binny_dir,"/binny_WS.Rdata"))
    print("DONE: Saving R image: binny_WS.Rdata")
  }}
