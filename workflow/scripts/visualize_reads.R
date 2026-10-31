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

###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input[[1]] #"Stats/all_stats.Rdata"
outplot             <- snakemake@output #"Visualization/{types}.reads.png"

###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

###################################################################################################
## Fields of interest
###################################################################################################
foi <- data.frame(
  "steps"=c("Preprocessing/mg.r1.fq","Preprocessing/mg.r1.preprocessed.fq",
            "Preprocessing/mg.se.preprocessed.fq",
            "Preprocessing/mg.r1.trimmed.fq",
            "Preprocessing/mg.se.trimmed.fq","Preprocessing/mt.r1.fq",
            "Preprocessing/mt.r1.preprocessed.fq",
            "Preprocessing/mt.se.preprocessed.fq","Preprocessing/mt.r1.trimmed.fq",
            "Preprocessing/mt.se.trimmed.fq",
            "Preprocessing/mt.r1.trimmed__rrna_proc__fq",
            "Preprocessing/mt.se.trimmed__rrna_proc__fq","Assembly/mg.reads.sorted.bam",
            "Assembly/mt.reads.sorted.bam",
            "kraken classified reads","mOTU mapping reads",
            "paired reads","se reads",
            "primary lines in bam","segments mapped","both partners map",
            "one partner maps",c(paste("mapping to",dbs))),
  "rank"=c(1,5:6,3:4,1,7:8,3:4,5:6,9,9,101,102,11:16,20+1:length(dbs)),
  "part"=c(rep("preprocessing",12),rep("assembly",2),rep("taxonomy",2),rep("preprocessing",2),
           rep("assembly",4),rep("annotation",length(dbs))),
  "name"=c("raw","preprocessed pairs","preprocessed singletons",
           "trimmed pairs","trimmed singletons","raw",
           "preprocessed pairs","preprocessed singletons",
           "trimmed pairs", "trimmed singletons",
           "filtered pairs","filtered singletons","reads in alignment",
           "reads in alignment",
           "kraken classified reads","mOTU mapping reads",
           "pairs in alignment","singletons in alignment",
           "primary lines in bam","segments mapped","both partners map",
           "one partner maps",c(paste("mapping to",dbs))),
  "mg_col"=c(colorRampPalette(brewer.pal(9,"Reds")[4:7])(6)[c(1,4,4,2,2,1,4,4,2,2,3,3,5,5)],
             brewer.pal(9,"Set1")[8],brewer.pal(8,"Set2")[4],
             colorRampPalette(brewer.pal(9,"Reds")[4:7])(6)[c(6,6)],
             brewer.pal(9,"Purples")[5:7][c(1,2,3,3)],
             colorRampPalette(brewer.pal(9,"Set1")[c(6,3)])(length(dbs))[1:length(dbs)]),
  "mt_col"=c(colorRampPalette(brewer.pal(9,"Blues")[4:7])(6)[c(1,4,4,2,2,1,4,4,2,2,3,3,5,5)],
             brewer.pal(9,"Set1")[8],brewer.pal(8,"Set2")[4],
             colorRampPalette(brewer.pal(9,"Blues")[4:7])(6)[c(6,6)],
             brewer.pal(9,"Purples")[5:7][c(1,2,3,3)],
             colorRampPalette(brewer.pal(9,"Set1")[c(6,3)])(length(dbs))[1:length(dbs)]),
  stringsAsFactors = F)

###################################################################################################
## plot
###################################################################################################
print("Plotting reads stats.")
dir.create("Visualization", showWarnings = FALSE)
for(type in types){
  readM <- merge(statsTab[statsTab$type==type&
                            statsTab$step %in% c(foi$steps),],
                 foi,all.x=T,by.x="step",by.y="steps")
  readM <- readM[order(readM$rank),]
  plotfile <- grep(paste0("/",type,"."),fixed=T,unlist(outplot),value = T)
  png(plotfile,width=8.5,height=7,units = "cm",res=300,pointsize=7)
  par("mar"=c(10,5,0.5,0.5),tcl="-0.3",mgp=c(3.5,0.5,0),lwd=0.5)
  a <- barplot(readM$number,
               las=2,names.arg = readM$name,
               ylim=c(0,max(readM$number)*1.09),yaxs="i",ylab="numbers",axes=F,
               col=sapply(type, function(x) if(x=="mg") readM$mg_col else readM$mt_col),
               border=NA,cex.lab=9/7)
  axis(2,lwd=0.5,las=1)
  abline(v=unlist(sapply(1:(nrow(readM)-1),function(x) if(readM$part[x] != readM$part[x+1]) (a[x]+a[x+1])/2)))
  box(bty="o")
  dev.off()
}
