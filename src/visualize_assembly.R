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

stats_image         <- snakemake@input #"Stats/all_stats.Rdata"
outplot             <- snakemake@output # "Visualization/Assembly_length_hist.png"
#                                         Visualization/Assembly_cov_hist.png"
                                        #"Visualization/Assembly_gc_hist.png"
                                         #"Visualization/Assembly_cov_gc_scatter.png"

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
            "Preprocessing/mt.r1.trimmed.rna_filtered.fq",
            "Preprocessing/mt.se.trimmed.rna_filtered.fq","Assembly/mg.reads.sorted.bam",
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
print("Plotting assembly stats.")
dir.create("Visualization", showWarnings = FALSE)
plotfile <- outplot[[1]]
length_hist_log <- hist(log10(contig_length$length),breaks=100,plot=F)
png(plotfile,width=5.7,height=4.6,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
plot(length_hist_log$mids,length_hist_log$counts,type="l",
          las=1,xlim=c(2,max(length_hist_log$mids)*1.03),xaxs="i",
          ylim=c(0,max(length_hist_log$counts)*1.04),yaxs="i",
          ylab="",col= brewer.pal(9,"Purples")[7],
          xlab=expression("log"["2"]*" contig length"),cex.lab=9/7,axes=F,lwd=0.8)
abline(v=log10(statsTab$number[statsTab$step=="N50"]),lty=3)
text(log10(statsTab$number[statsTab$step=="N50"]),
     max(length_hist_log$counts)*0.95,
     "labels"=paste("N50 =",statsTab$number[statsTab$step=="N50"]),pos = 4)
mtext("counts",2,2.8,cex=9/7)
axis(1,at = 2:6,lwd=0.5,las=1)
axis(2,lwd=0.5,las=1)
box(bty="o")
dev.off()

plotfile <- outplot[[2]]
png(plotfile,width=5.7,height=4.6,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
depth_hist_log <- list()
for(type in types){
  depth_hist_log[[type]] <- hist(log10(depthAssStat[[type]]$depth),
                         breaks=seq(from=log10(min(sapply(depthAssStat,function(x) min(x$depth)))),
                                    to=log10(max(sapply(depthAssStat,function(x) max(x$depth)))),
                                    length.out = 100),plot=F)
}
plot(1,type="n",axes=F,ylab="",xlab=expression("log"["2"]*" contig coverage depth"),cex.lab=9/7,
     xlim=c(min(depth_hist_log[[1]]$mids)*ifelse(min(depth_hist_log[[1]]$mids)<0,1.03,0.97),
            max(depth_hist_log[[1]]$mids)*1.03),xaxs="i",
     ylim=c(0,max(sapply(depth_hist_log,function(x) max(x$counts)))))
mtext("counts",2,2.8,cex=9/7)
axis(1,lwd=0.5,las=1)
axis(2,lwd=0.5,las=1)
box(bty="o")
for(type in types){
  lines(depth_hist_log[[type]]$mids,depth_hist_log[[type]]$counts,type="l",
        col =ifelse(type=="mg",brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7]),
        lwd=0.8)
}
dev.off()

plotfile <- outplot[[3]]
gc_hist <- hist(contig_gc$GC,breaks=100,plot=F)
png(plotfile,width=5.7,height=4.6,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.65),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
plot(gc_hist$mids,gc_hist$counts,type="l",
     las=1,xlim=c(0,100),xaxs="i",
     ylim=c(0,max(gc_hist$counts)*1.04),yaxs="i",
     ylab="",col= brewer.pal(9,"Purples")[7],
     xlab="contig %GC",cex.lab=9/7,axes=F,lwd=0.8)
mtext("counts",2,2.8,cex=9/7)
axis(1,lwd=0.5,las=1)
axis(2,lwd=0.5,las=1)
box(bty="o")
dev.off()

plotfile <- outplot[[4]]
if(length(types)>1){
  spec <- colorRampPalette(rainbow(9)[9:1])
  plotcol <- spec(length(unique(perContM$GC)))[as.numeric(as.factor(perContM$GC))]
  
  perContM <-merge(depthAssStat[[1]],depthAssStat[[2]],by=1,all=T,suffixes=paste0(".",types))
  perContM <-merge(perContM,contig_gc,by=1,all=T)
  perContM[is.na(perContM)] <- 0
  png(plotfile,width=8.5,height=8.5,units = "cm",res=300,pointsize=7)
  par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  plot(perContM$depth.mg,perContM$depth.mt,
       pch=16,cex=0.2,log="xy",
       las=1,xlim=c(min(depthAssStat[["mg"]]$depth)/11,max(depthAssStat[["mg"]]$depth)/1.2),xaxs="i",
       ylim=c(min(depthAssStat[["mt"]]$depth)/11,max(depthAssStat[["mt"]]$depth)/1.2),yaxs="i",
       ylab="",
       xlab="contig coverage depth (metaG)",cex.lab=9/7,axes=F,
       col=plotcol)
  points(perContM$depth.mg[perContM$depth.mt==0],
         rep(min(depthAssStat[["mt"]]$depth)/10,length(which(perContM$depth.mt==0))),
         pch=16,cex=0.2,
         col=plotcol[perContM$depth.mt==0])
  points(rep(min(depthAssStat[["mg"]]$depth)/10,length(which(perContM$depth.mg==0))),
         perContM$depth.mt[perContM$depth.mg==0],
         pch=16,cex=0.2,
         col=plotcol[perContM$depth.mg==0])
  mtext("contig coverage depth (metaT)",2,2.8,cex=9/7)
  axis(1,
       at=c(min(depthAssStat[["mg"]]$depth)/10,10^c(-3:8)),
       labels=c(0,format(10^c(-3:8),scientific=T,digits = 2)),lwd=0.5,las=1)
  axis(2,at=c(min(depthAssStat[["mt"]]$depth)/10,10^c(-3:8)),
       labels=c(0,format(10^c(-3:8),scientific=T,digits = 2)),lwd=0.5,las=1)
  box(bty="o")
  dev.off()
}else{
  perContM <-merge(depthAssStat[[1]],contig_gc,by=1,all=T)
  perContM[is.na(perContM)] <- 0
  png(plotfile,width=8.5,height=8.5,units = "cm",res=300,pointsize=7)
  par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  plot(perContM$GC,perContM$depth,
       pch=16,cex=0.2,log="y",
       las=1,ylim=c(min(depthAssStat[[1]]$depth)/11,max(depthAssStat[[1]]$depth)/1.2),yaxs="i",
       xlim=c(0,100),xaxs="i",
       ylab="",
       xlab="contig %GC",cex.lab=9/7,axes=F,
       col=brewer.pal(9,"Purples")[7])
  if(any(perContM$depth==0)){
    points(perContM$GC[perContM$depth==0],
           rep(min(depthAssStat[[1]]$depth)/10,length(which(perContM$depth==0))),
           pch=16,cex=0.2,
           col=brewer.pal(9,"Purples")[7])
  }
  mtext(paste0("contig coverage depth (meta",ifelse(types[1]=="mg","G","T"),")"),
        2,2.8,cex=9/7)
  axis(1,lwd=0.5,las=1)
  axis(2,lwd=0.5,las=1)
  box(bty="o")
  dev.off()
}
