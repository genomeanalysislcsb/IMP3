#!/bin/R

###################################################################################################
## Load required packages
###################################################################################################

print("Loading required R libraries")
require(genomeIntervals)

#require(checkpoint)
###################################################################################################
## Read in arguments
###################################################################################################
# User should input the three different files
# Uncomment when ready

print("START: Reading arguments")

args                <- commandArgs(trailingOnly = TRUE)
locations	    <- args[1]
out_dir             <- args[2] # Output directory
MG.read.count_file  <- args[3]
MG.map.summary_file <- args[4]
MG.cov_file         <- args[5]
MG.depth_file       <- args[6]
#MG.var_file         <- args[7]
GC.dat_file         <- args[7]
coords_file         <- args[8]
annot_file          <- args[9]
#nucmer_file         <- args[11]
function_script     <- args[10]
pk 		    <- as.numeric(args[11])
nn 		    <- as.numeric(args[12])

print("DONE: Reading arguments")

#checkpoint('2016-06-20', scanForPackages=FALSE, checkpointLocation=locations, project=locations)

require(ggplot2)
require(gtools)
require(gtable)
require(data.table)
require(reshape)
require(grid)
require(grDevices)
require(stringr)
require(xtable)
require(beanplot)
require(psych)

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
## Initialize functions for various calculations and normalizations
###################################################################################################
print("START: Reading functions")
source(function_script)
print("DONE: Reading functions")

###################################################################################################
## Read in the necessary input files
###################################################################################################
print("START: Reading data")

## Read counts MG
print("Read in MG read count file")
MG.read.count <- read.table(MG.read.count_file)
MG.read.count <- as.data.frame(cbind(as.character(MG.read.count$V2)[-nrow(MG.read.count)], MG.read.count$V1[-nrow(MG.read.count)]))
colnames(MG.read.count) <- c("file", "count")
MG.read.count$file <- as.character(MG.read.count$file)
# Divide by four for total sequences because original file is a line count
MG.read.count$count <- as.numeric(as.character(MG.read.count$count))/4

## Read MG mapping summary file (samtools flagstat)
print("Read in MG mapping summary file")
MG.map.summary <- read.table(MG.map.summary_file)
MG.map.summary <- cbind(
		     c("Total",
		       "Duplicates",
		       "Mapped",
		       "Paired",
		       "Read1",
		       "Read2",
		       "Properly paired",
		       "With itself and mate",
		       "Singletons",
		       "Mate mapped to different contig",
		       "Mate mapped to different contig mapQ>=5"
		       ),
		     MG.map.summary
		     )
colnames(MG.map.summary) <- c("Mapping", "Reads")


## MG coverage (samtools idxstat output)
print("Read in MG coverage file")
MG.cov <- read.table(MG.cov_file, colClasses=c("factor", rep("numeric", 6 )))[,c(1,4:7)]
colnames(MG.cov) <- c("contig", "MG_reads", "MG_bases_covered", "length", "MG_cov")

## MG depth (bedtools coverageBed output)
print("Read in MG depth file")
MG.depth <- read.table(MG.depth_file, colClasses=c("factor", "numeric"))
colnames(MG.depth) <- c("contig", "MG_depth")

## MG variant calls
# here should be a filterstep - AHB
#print("Read in MG variant calls")
#MG.var <- as.data.frame(table(read.table(MG.var_file)[,1]))
#colnames(MG.var) <- c("contig", "MG_var")

## GC content file
print("Read in GC percentage file")
GC.dat <- read.delim(GC.dat_file)
colnames(GC.dat) <- c("contig", "GC")

## gff annotation file
print("Read in gff3 annotation file")
annot <- readGff3(annot_file, isRightOpen = T)

# extract only interesting information
# this needs to be adapted to the new format -AHB
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
#annot.1$database[annot.1$database==""] <- NA
annot.1[,2:3] <- sapply(annot.1[2:3], as.character)
annot.1[,2:3] <- sapply(annot.1[2:3], as.numeric)



# calculate no. of gene length
print("Calculating no. of genes per contig")
annot.2 <- as.data.frame(table(annot.1$contig))
colnames(annot.2)[1:ncol(annot.2)] <- c("contig","no_of_genes")

# calculate gene lengths
print("Calculating length of genes from gff3 file")
annot.1 <- as.data.frame(cbind(annot.1, annot.1$end - annot.1$start + 1))
colnames(annot.1)[ncol(annot.1)] <- "gene_length"

# aggregate table and calculate total gene lengths within contig
save.image(name_plot("mg_results.Rdat"))
print("Calculating coding density of contigs")
# Create temporary table
total_gene_length <- aggregate(gene_length~contig, data=annot.1, FUN=sum)
colnames(total_gene_length)[ncol(total_gene_length)] <- "total_gene_length"
annot.2 <- merge(annot.2,
		 total_gene_length,
		 by="contig",
		 all.y=FALSE)

# remove temporary table
rm(total_gene_length)

# create annotation table
print("Creating annotation table")
annot.E.prelim <- annot.1[,c("contig","start","end","essential")]
annot.E.prelim <- annot.E.prelim[!is.na(annot.E.prelim$essential),]
annot.E.prelim$essential <- as.character(annot.E.prelim$essential)
annot.E.prelim$contig <- as.character(annot.E.prelim$contig)
for(column in dbs) annot.1[,column] <- as.character(annot.1[,column])
annot.1[,dbs][is.na(annot.1[,dbs])]  <- 0
annot.1[,dbs][annot.1[,dbs]!=0]  <- 1
for(column in dbs) annot.1[,column] <- as.numeric(annot.1[,column])

annot.3 <- aggregate(annot.1[,dbs],list(annot.1$contig),sum)
colnames(annot.3)[1] <- "contig"
annot.3$all_annotations <- rowSums(annot.3[,dbs])

annot_final <- merge(annot.2, annot.3, by="contig")

# vizbin points, the contig names are usually not included in this file,
# concatenate it before reading in
print("Reading in vizbin coordinates")
coords <- read.table(coords_file, colClasses=c("factor", "numeric", "numeric"),
		     sep="\t", col.names=c("contig", "x", "y"))

# Read in nucmer output from metaquast
#print("Reading in nucmer results")
#nucmer_try <- try(read.table(nucmer_file, header=F), silent=T)
#if(inherits(nucmer_try, "try-error")){
#  print("WARNING: Nucmer file empty. No taxonomy was assigned to contigs")
#  nucmer_res <- read.table(text = "",
#			   col.names = c("ref_start", "ref_end", "query_start", "query_end",
#					 "ref_align_len", "query_align_len", "identity",
#					 "ref_id", "contig")
#	     )
#}else{
#  nucmer_res <- read.table(nucmer_file, header=F)
#  colnames(nucmer_res) <- c("ref_start", "ref_end", "query_start", "query_end", "ref_align_len",
#			  "query_align_len", "identity", "ref_id", "contig")
#}

print("DONE: Reading data")
###################################################################################################
## Merge the data sets without the vizbin coordinates
###################################################################################################
print("START: Merging data")

all.dat <- merge(MG.cov, GC.dat, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, MG.depth, by=c("contig"), all=T, incomparables=NA)
#all.dat <- merge(all.dat, MG.var, by=c("contig"), all=T, incomparables=NA)
all.dat <- merge(all.dat, annot_final, by=c("contig"), all=T, incomparables=NA)
#all.dat <- merge(all.dat, nucmer_res, by=c("contig"), all=T, incomparables=NA)

print("DONE: Merging data")
###################################################################################################
## Perform calculations and append it to the full table
###################################################################################################
# Calculate and merge data
print("START: Performing calculations")

# Compute gene density
print("Computing gene density")
all.dat$gene_dens <- gene_density(all.dat$no_of_genes, all.dat$length)

# Compute coding density
print("Computing coding density")
all.dat$coding_dens <- coding_density(all.dat$total_gene_length, all.dat$length)

# Compute RPKM values
print("Computing MG RPKM values")
all.dat$MG_rpkm <- contig_rpkm(all.dat$MG_reads, all.dat$length)

# Compute variant densities
#print("Computing MG variant densities")
#all.dat$MG_var_dens <-  var_density(all.dat$MG_var, all.dat$length, all.dat$MG_reads)

# Compute query coverage against reference genomes
#print("Computing query coverage of contigs against reference genomes")
#all.dat$query_cov <- all.dat$query_align_len / all.dat$length * 100
#all.dat$query_cov[na.omit(all.dat$query_cov > 100)] = 100

save.image(name_plot("mg_results.Rdat"))
print("DONE: Performing calculations")
###################################################################################################
## Organize filtering statistics and create table
###################################################################################################
print("START: Organizing preprocessing stats")

## Organize MG filtering procedures
print("Printing metagenomic filtering statistics")
# Obtain file names without path
MG.read.count$file <- unlist(lapply(MG.read.count$file, get_file_name))

# Obtain string length
MG.read.count <- cbind(MG.read.count,
		       unlist(lapply(MG.read.count$file, nchar)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "string_len"

# Order the file names according to the string length
MG.read.count <- MG.read.count[order(MG.read.count$string_len),]

# Create additional column for fastq file types (paired or singletons)
MG.read.count <- cbind(MG.read.count, rep(NA, nrow(MG.read.count)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "type"
MG.read.count$type[grep("R[12]", MG.read.count$file)] <- "paired-end"
MG.read.count$type[grep("SE", MG.read.count$file)] <- "singletons"

# Create new columns with filtering type (within file names)
MG.read.count <- cbind(MG.read.count,
		      unlist(lapply(MG.read.count$file, filtering)))
colnames(MG.read.count)[ncol(MG.read.count)] <- "filtering"
MG.read.count$filtering <- as.character(MG.read.count$filtering)

# Rename filtering steps
MG.read.count$filtering[ # Unfiltered/raw
			grep("R[12]", MG.read.count$filtering)] <-
			   "unfiltered (raw)"
MG.read.count$filtering <- # Remove _filt extension in file name
   gsub("_filt", "", MG.read.count$filtering, ignore.case=TRUE)
MG.read.count$filtering <- # Remove _filt extension in file name
   gsub("trimmed", "adapter trimmed and quality filtered", MG.read.count$filtering, ignore.case=TRUE)
MG.read.count$filtering <- # Convert uniq to unique
   gsub("uniq", "collapsed unique", MG.read.count$filtering, ignore.case=TRUE)
#MG.read.count$filtering <- # Convert hg to "human sequences"
#    gsub("hg\\d+", "human sequences filtered", MG.read.count$filtering,
#	ignore.case=TRUE)

# Summarize table for final report
MG.read.count.final <- unique(MG.read.count[,c(5,4,2)])
MG.read.count.final$count <- as.integer(MG.read.count.final$count)

# Write out table in html and tab separated files
sink(name_plot("mg.read_stats.html"))
print(xtable(MG.read.count.final, html.table.attributes=""), type = "html")
sink()

write.table(MG.read.count.final, name_plot("mg.read_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

print("DONE: Organizing preprocessing stats")

###################################################################################################
## Incorporate coordinates from vizbin
###################################################################################################
print("START: Incorporating VizBin data")

vb_dat <- merge(all.dat, coords, by=c("contig"), all=F, incomparables=NA)
vb_dat <- vb_dat[!is.na(vb_dat$x),]

print("Processing VizBin data")
# Handle missing and infinite values
vb_dat[is.na(vb_dat)] <- 0


#binny
LIB <- "Binning"
contigInfo <- vb_dat
contigInfo$contig <- as.character(contigInfo$contig)
annot.E <- annot.E.prelim[annot.E.prelim$contig %in% vb_dat$contig,]
essAnno <- aggregate(annot.E$essential,list(annot.E$contig),function(x) paste(x,sep=";",collapse=";"))
colnames(essAnno) <- c("contig","essentialGene")
contigInfo <- merge(contigInfo,essAnno,by.x=1,by.y=1,all.x=T)
contigInfo$essentialGene[is.na(contigInfo$essentialGene)] <- "notEssential"

# auto clustering
coco <- c(which(colnames(contigInfo)=="x"),which(colnames(contigInfo)=="y"))
#estimate reachability distance (based on number of neighbouting points)
skn <- sort(knn.dist(contigInfo[,coco],pk)[,pk])
sdkn <- runsd(skn,10)
est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))]
#document estimated reachability distances
write.table(t(c("scan","reachabilityDist")),paste("Binning/binny/reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t")
write.table(t(c("first",est)),paste("Binning/binny/reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
#DBscan
cdb <- dbscan(contigInfo[,coco],est,pk)

#count contigs in clusters and document clustering steps
pdf(paste("Binning/binny/scatterPlot1",pk,pk,"pdf",sep="."))
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

write.table(cdbTab,paste("Binning/binny/clusterFirstScan",pk,nn,"tsv", sep="."),sep="\t",row.names=F,quote=F)
write.table(t(c("clusterName","cutoff")),paste("Binning/binny/bimodalClusterCutoffs",pk,nn,"tsv",sep="."),sep="\t",row.names=F,col.names=F,quote=F)

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
write.table(t(c("cluster","cutoff")),paste("Binning/binny/bimodalClusterCutoffs",pk,nn,"tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
for(ds in unique(clusterRes$cluster[grep("D",clusterRes$cluster)])){
  clusterRes <- muClus(ds,clusterRes,paste("Binning/binny/bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
}

#second iteration of clustering on clusters with more than 20% duplicated essential genes 
### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
### - neighbouring points are increased by 2
pk2 <- pk + 2
pdf(paste("Binning/binny/scatterPlots2",pk,nn,"pdf",sep="."))
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

      write.table(t(c(bb,estBB)),paste("Binning/binny/reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
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
      write.table(BBcdbTab,paste("Binning/binny/clusterFiles/blob1Cluster",pk,nn,bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)

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
        clusterRes <- muClus(bds,clusterRes,paste("Binning/binny/bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
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
pdf(paste("Binning/binny/scatterPlots3",pk,nn,"pdf",sep="."))
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

      write.table(t(c(bb,estBB)),paste("Binning/binny/reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
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
      write.table(BBcdbTab,paste("Binning/binny/clusterFiles/blob2Cluster",pk,nn,bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)

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
        clusterRes <- muClus(bds,clusterRes,paste("bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
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
pdf(paste("Binning/binny/scatterPlots4",pk,nn,"pdf",sep="."))
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

      write.table(t(c(bb,estBB)),paste("Binning/binny/reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
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
      write.table(BBcdbTab,paste("Binning/binny/clusterFiles/blob3Cluster",pk,nn,bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)

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
        clusterRes <- muClus(bds,clusterRes,paste("Binning/binny/bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
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
write.table(clusterRes,paste("Binning/binny/contigs2clusters",pk,nn,"tsv",sep="."),sep="\t",row.names=F,quote=F)
save(list=c("contigInfo","coco","pk","nn","essAnno","LIB","find.cutoff","muClus","clusterRes"),file=paste("Binning/binny/clusteringWS",pk,nn,"Rdata",sep="."))
saveRDS(clusterRes,paste("Binning/binny/contigs2clusters",pk,nn,"RDS",sep="."))

print("Saved tables. Plotting final map now.")
#plot final map
essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(111)[111:1]
pdf(paste("Binning/binny/finalClusterMap",pk,nn,"pdf",sep="."),width=6,height=5,pointsize=8)
layout(matrix(1:2,nrow=1,ncol=2),c(5/6,1/6))
par(mar=rep(1,4),pty="s")
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F,bty="o")
cI <- data.frame("essNum"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "medX"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "medY"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "high"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
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
mtext(LIB,3,0)
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

rm(list=c("contigInfo","coco","pk","nn","LIB","clusterRes","essPal","cI","cI2","hi","a","clus","ess","BBduClus","BBclusterRes","BBduCont","uniD","BBcloseClus","BBcloseCont","BBemptyClus","BBemptyCont","BBcdbTab","bbInfo","estBB","BBcdb","sknBB","emptyClus","closeClus","cdbTab","skn","sdkn","est","sdknBB","pk2","pk4","pk6"))


# Remove outliers and infinite values
print("Handling outliers")

#vb_dat$MG_var_dens <- outliers(vb_dat$MG_var_dens,2)
vb_dat$MG_depth <- outliers(vb_dat$MG_depth,2)
vb_dat$MG_rpkm <- outliers(vb_dat$MG_rpkm,2)

save.image(name_plot("mg_results.Rdat"))
print("DONE: Incorporating VizBin data")


####################################################################
## ASSEMBLY STATISTICS AND VISUALIZATIONS
####################################################################
print("Begin visualizing data")
## Output filtering statistics table (text and tsv)

## Output assembly statistics table (text and tsv)
print("Print assembly statistics table")
sink(name_plot("assembly_stats.html"))
print(xtable(assembly.stats, html.table.attributes=""), type = "html")
sink()

write.table(assembly.stats, name_plot("assembly_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

## Output mapping stats table
print("Print metagenomic mapping statistics table")
sink(name_plot("mg_mapping_stats.html"))
print(xtable(MG.map.summary, html.table.attributes=""), type="html")
sink()

write.table(MG.map.summary, name_plot("mg_mapping_stats.txt"),
	    sep="\t", quote=F,
	    row.names=F)

## Plot standard vizbin scatter plot (non included)
print("Generating standard vizbin plot")
png(name_plot("IMP-vizbin_standard.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point() +
ggtitle("Standard VizBin") +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length info
print("Generating vizbin plot with contig length information")
png(name_plot("IMP-vizbin_length.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length))) +
guides(size=guide_legend(title=log10len)) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and GC content
print("Generating vizbin plot for GC content")
png(name_plot("IMP-vizbin_length_GC.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(aes(size=log10(length), colour=GC, order=GC),alpha=0.5) +
scale_colour_gradientn(colours=rainbow(7), guide="colourbar", guide_legend(title="% G+C")) +
guides(size=guide_legend(title=log10len)) +
theme_nothing()
dev.off()

print("OVER")
####################################################################
## MAPPING STATISTICS AND VISUALIZATIONS
####################################################################

## Plot mappable reads density
print("Generating mapped reads plot")
var1 <-log10(c(all.dat$MG_reads,all.dat$MG_rpkm))
var1[is.infinite(var1)]=NA
var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
MG_mapped_reads<-data.frame(var1,var2)

png(name_plot("IMP-mg_reads_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MG_mapped_reads,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Mappable reads", ylab=expression(log[10]*~"count"))
legend("bottomleft", fill =c("blue", "red"), legend = c("No of reads mapped", "RPKM normalized"))
dev.off()

## Plot coverage density
print("Generating MG coverage plot")
var1 <-c(all.dat$MG_cov,all.dat$MG_depth)
var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
MG_coverage<-data.frame(var1,var2)

png(name_plot("IMP-mg_coverage_density.png"), width=350, height=700)
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
beanplot(var1 ~ var2, data= MG_coverage,  side = "both",log="auto",
what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
bw="nrd0", main="Coverage", ylab="fraction")
legend("bottomleft", fill =c("blue", "red"), legend = c("Coverage", "depth"))
dev.off()

## Plot vizbin scatter plot with length and MG coverage info
print("Generating vizbin plot for metagenomic coverage")
png(name_plot("IMP-mg_vizbin_length_cov.png"), width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_cov, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "coverage"))))
      ) +
theme_nothing()
dev.off()

## Plot vizbin scatter plot with length and MG depth info
print("Generating vizbin plot for metagenomic depth")
png(name_plot("IMP-mg_vizbin_length_depth.png"),width=700, height=700)
ggplot(vb_dat, aes(x=x,y=y)) +
geom_point(colour="blue", aes(alpha=MG_depth, size=log10(length))) +
guides(size=guide_legend(title=log10len),
       alpha=guide_legend(title=expression(bold(atop("Metagenomic", "depth"))))
      ) +
theme_nothing()
dev.off()

####################################################################
## VARIANT STATISTICS AND VISUALIZATIONS
####################################################################
## Plot variant count
#var1 <-log10(c(all.dat$MG_var,all.dat$MG_dens))
#var1[is.infinite(var1)]=NA
#var2 <- c(rep("MG",nrow(all.dat)),rep("MG",nrow(all.dat)))
#MG_variant_count<-data.frame(var1,var2)

#print("Generating variant count plots")
#png(name_plot("IMP-mg_var_count.png") ,width=350, height=700)

#par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
#beanplot(var1 ~ var2, data= MG_variant_count,  side = "both",log="auto",
#what=c(1,1,1,0), border = NA, col = list("blue", c("red", "white")),
#bw="nrd0", main="MG variant (SNPs & INDELS)", ylab=expression(log[10]*~count))
#legend("bottomleft", fill =c("blue", "red"), legend = c("No. of variants", "Variant density"))
#dev.off()

## Plot vizbin scatter plot with length and MG variant info
# Create label
#MG_var_label <- expression(bold(frac(variants[MG]/kb, MG[rpkm])))

#print("Generating vizbin plot for metagenomic variant density")
#png(name_plot("IMP-mg_vizbin_length_vardens.png"), width=700, height=700)
#ggplot(vb_dat, aes(x=x,y=y)) +
#geom_point(aes(colour=MG_var_dens, size=log10(length), order=MG_var_dens), alpha=0.75) +
#scale_colour_gradient(high="black", low="cyan") +
#       guides(size=guide_legend(title=log10len),
#	      colour=guide_colourbar(title=MG_var_label)
#       ) +
#theme_nothing()
#dev.off()


####################################################################
## Save the R workspace
####################################################################

print("START: Saving R image: mg_results.Rdat")
save.image(name_plot("mg_results.Rdat"))
print("DONE: Saving R image: mg_results.Rdat")



