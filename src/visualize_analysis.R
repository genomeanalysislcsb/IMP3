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
library(UpSetR)

###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input[[1]] #"Stats/all_stats.Rdata"
outplot             <- snakemake@output # Visualization/Analysis_length_hist.png"
                                        #"Visualization/Analysis_cov_hist.png
                                        # Visualization/Analysis_cov_scatter.png
                                        # Visualization/Analysis_func_sets.png

###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

expand_hist_mids <- function(hist){
  binsize <- abs(hist$mids[2]-hist$mids[1])
  mids <- c(hist$mids-binsize,hist$mids[length(hist$mids)],hist$mids[length(hist$mids)]+binsize)
  return(mids)
}

expand_hist_counts <- function(hist){
  binsize <- abs(hist$mids[2]-hist$mids[1])
  counts <- c(0,hist$counts,0)
  return(counts)
}

gff2complete <- function(gff){
  completeness <- data.frame(t(sapply(paste0(gff$attribute[gff$feature=="CDS"],";"), 
                                      function(x) c(gsub("ID=([^;]+);.+","\\1",x),
                                                    gsub(".+partial=([01]{2});.*","\\1",x)))),
                                                        stringsAsFactors = F)
  colnames(completeness) <- c("gene","partial")
  return(completeness)                        
}

gff2func <- function(gff,func){
  funcT <- data.frame(t(sapply(paste0(gff$attribute[gff$feature=="CDS"&
                                               grepl(paste0(".+;",func,"=([^;]+)"),gff$attribute)],";"), 
                                      function(x) c(gsub("ID=([^;]+);.+","\\1",x),
                                                    gsub(paste0(".+;",func,"=([^;]+);.*"),"\\1",x)))),
                             stringsAsFactors = F)
  colnames(funcT) <- c("gene",func)
  return(funcT)                        
}

min0 <- function(x){
  min(x[x>0])
}

###################################################################################################
## plot
###################################################################################################
completeness <- gff2complete(gene_gff)

print("Plotting analysis stats.")
dir.create("Visualization", showWarnings = FALSE)
plotfile <- outplot[[1]]
length_hist_log <- hist(log10(gene_length$length),breaks=100,plot=F)
length_comp_hist_log <- hist(log10(gene_length$length[gene_length$gene %in% completeness$gene[completeness$partial=="00"]]),
                             breaks=length_hist_log$breaks,plot=F)
length_part_hist_log <- hist(log10(gene_length$length[gene_length$gene %in% completeness$gene[completeness$partial!="00"]]),
                             breaks=length_hist_log$breaks,plot=F)

png(plotfile,width=5.7,height=4.6,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
plot(expand_hist_mids(length_hist_log),
     expand_hist_counts(length_hist_log),type="l",
          las=1,xlim=c(log10(50),max(length_hist_log$mids)*1.03),xaxs="i",
          ylim=c(0,max(length_hist_log$counts)*1.04),yaxs="i",
          ylab="",col= brewer.pal(9,"RdPu")[7],
          xlab=expression("log"["2"]*" gene length"),cex.lab=9/7,axes=F,lwd=0.8)
lines(expand_hist_mids(length_comp_hist_log),
     expand_hist_counts(length_comp_hist_log),col=brewer.pal(9,"RdPu")[7],lwd=1)
lines(expand_hist_mids(length_part_hist_log),
      expand_hist_counts(length_part_hist_log),col=brewer.pal(9,"RdPu")[7],lwd=0.8,lty=3)
mtext("counts",2,2.8,cex=9/7)
axis(1,at = 0:6,lwd=0.5,las=1)
axis(2,lwd=0.5,las=1)
box(bty="o")
dev.off()

plotfile <- outplot[[2]]
png(plotfile,width=5.7,height=4.6,units = "cm",res=300,pointsize=7)
par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
depth_hist_log <- list()
for(type in types){
  depth_hist_log[[type]] <- hist(log10(depthGeneStat[[type]]$depth),
                         breaks=seq(from=log10(min(sapply(depthGeneStat,function(x) min0(x$depth)))),
                                    to=log10(max(sapply(depthGeneStat,function(x) max(x$depth)))),
                                    length.out = 100),plot=F)
}
plot(1,type="n",axes=F,ylab="",xlab=expression("gene coverage depth"),cex.lab=9/7,
     xlim=c(expand_hist_mids(depth_hist_log[[1]])[1]-1.03,
            max(expand_hist_mids(depth_hist_log[[1]]))*1.03),xaxs="i",
     ylim=c(0,max(c(sapply(depth_hist_log,function(x) max(x$counts)),
                  sapply(depthGeneStat,function(x) length(which(x$depth==0)))))*1.1),yaxs="i")
mtext("counts",2,2.8,cex=9/7)
axis(1,at=c(expand_hist_mids(depth_hist_log[[1]])[1]-1,-3:10),
     labels = format(c(0,10^(-3:10)),scientific=T),lwd=0.5,las=1)
axis(2,lwd=0.5,las=1)
for(type in types){
  lines(expand_hist_mids(depth_hist_log[[type]]),
        expand_hist_counts(depth_hist_log[[type]]),type="l",
        col =ifelse(type=="mg",brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7]),
        lwd=0.8)
  points(expand_hist_mids(depth_hist_log[[1]])[1]-ifelse(type=="mg",1,0.98),
         length(which(depthGeneStat[[type]]$depth==0)),
         type="h",col=ifelse(type=="mg",brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7]))
  text(expand_hist_mids(depth_hist_log[[1]])[1]-ifelse(type=="mg",1,0.98),
       length(which(depthGeneStat[[type]]$depth==0)),
       labels = paste("\n",length(which(depthGeneStat[[type]]$depth==0)),"genes\nwithout coverage"),
       col=ifelse(type=="mg",brewer.pal(9,"Reds")[7],brewer.pal(9,"Blues")[7]),
       pos=4,cex=0.6)
}
box(bty="o")
dev.off()

func <- lapply(dbs, function(x) gff2func(gene_gff,x))

plotfile <- outplot[[3]]
if(length(types)>1){
  perGeneM <-merge(depthGeneStat[[1]],depthGeneStat[[2]],by=1,all=T,suffixes=paste0(".",types))
  perGeneM[is.na(perGeneM)] <- 0
  perGeneM <- merge(perGeneM, func[[which.max(sapply(func,nrow))]],
                         by=1,all.x=T)
  perGeneM[is.na(perGeneM)] <- ""
  plotcol <- c("grey",rainbow(length(unique(perGeneM[,4]))-1))[as.numeric(as.factor(perGeneM[,4]))]
  
  png(plotfile,width=8.5,height=8.5,units = "cm",res=300,pointsize=7)
  par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  plot(perGeneM$depth.mg,perGeneM$depth.mt,
       pch=16,cex=0.2,log="xy",
       las=1,xlim=c(min0(depthGeneStat[["mg"]]$depth)/11,max(depthGeneStat[["mg"]]$depth)/1.2),xaxs="i",
       ylim=c(min0(depthGeneStat[["mt"]]$depth)/11,max(depthGeneStat[["mt"]]$depth)/1.2),yaxs="i",
       ylab="",
       xlab="contig coverage depth (metaG)",cex.lab=9/7,axes=F,
       col=plotcol)
  points(perGeneM$depth.mg[perGeneM$depth.mt==0],
         rep(min0(depthGeneStat[["mt"]]$depth)/10,length(which(perGeneM$depth.mt==0))),
         pch=16,cex=0.2,
         col=plotcol[perGeneM$depth.mt==0])
  points(rep(min0(depthGeneStat[["mg"]]$depth)/10,length(which(perGeneM$depth.mg==0))),
         perGeneM$depth.mt[perGeneM$depth.mg==0],
         pch=16,cex=0.2,
         col=plotcol[perGeneM$depth.mg==0])
  mtext("gene coverage depth (metaT)",2,2.8,cex=9/7)
  axis(1,
       at=c(min0(depthGeneStat[["mg"]]$depth)/10,10^c(-3:8)),
       labels=c(0,format(10^c(-3:8),scientific=T,digits = 2)),lwd=0.5,las=1)
  axis(2,at=c(min0(depthGeneStat[["mt"]]$depth)/10,10^c(-3:8)),
       labels=c(0,format(10^c(-3:8),scientific=T,digits = 2)),lwd=0.5,las=1)
  box(bty="o")
  dev.off()
}else{
  perGeneM <-merge(depthGeneStat[[1]],func[[which.max(sapply(func,nrow))]],
                    by=1,all.x=T)
  perGeneM[is.na(perGeneM)] <- ""
  plotcol <- c("grey",rainbow(length(unique(perGeneM[,3]))-1))[as.numeric(as.factor(perGeneM[,3]))]
  
  png(plotfile,width=8.5,height=8.5,units = "cm",res=300,pointsize=7)
  par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  plot(1:nrow(perGeneM),
       sort(perGeneM$depth,decreasing = T),
       pch=16,cex=0.2,log="y",
       las=1,
       xlim=c(0,nrow(perGeneM)+1),xaxs="i",
       ylim=c(min0(depthGeneStat[[types[1]]]$depth)/11,max(depthGeneStat[[types[1]]]$depth)/1.2),yaxs="i",
       ylab="",
       xlab="rank",cex.lab=9/7,axes=F,
       col=plotcol[order(perGeneM$depth,decreasing = T)])
  points((length(which(perGeneM$depth>0))+1):nrow(perGeneM),
         rep(min0(depthGeneStat[[types[1]]]$depth)/10,length(which(perGeneM$depth==0))),
         pch=16,cex=0.2,
         col=plotcol[perGeneM$depth==0])
  mtext(paste0("contig coverage depth (meta",ifelse(types[1]=="mg","G","T"),")"),
        2,2.8,cex=9/7)
  axis(1,lwd=0.5,las=1)
  axis(2,at=c(min0(depthGeneStat[[1]]$depth)/10,10^c(-3:8)),
       labels=c(0,format(10^c(-3:8),scientific=T,digits = 2)),lwd=0.5,las=1)
  box(bty="o")
  dev.off()
}

plotfile <- outplot[[4]]
if(length(dbs)>1){
  func[["no_anno"]] <- data.frame("gene"=setdiff(gene_length$gene,unlist(sapply(func,function(x)x$gene))),
                                  "no_anno"="no_anno",stringsAsFactors = F)
  ff <- do.call(rbind,lapply(1:length(func), function(x) cbind(func[[x]]$gene,c(dbs,"no_anno")[x])))
  ulist <- tapply(rep(1,nrow(ff)),list(ff[,1],ff[,2]),function(y)y>0)
  ulist[is.na(ulist)] <-0
  png(plotfile,width=17,height=10,units = "cm",res=300,pointsize=7)
  par("mar"=c(2.9,4,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  upset(data.frame(ulist)[,rank(c(dbs,"no_anno"))],length(dbs)+1,100,
        sets.bar.color = c(colorRampPalette(brewer.pal(9,"Set1")[c(6,3)])(length(dbs)),"grey")[(1+length(dbs)):1],
        mainbar.y.label = "intersection gene number",sets.x.label = "annotated genes",keep.order = T,
        mb.ratio = c(0.5,0.5),order.by=c("degree","freq"))
  dev.off()
}else{
  png(plotfile,width=5.7,height=5.7,units = "cm",res=300,pointsize=7)
  par("mar"=c(0.5,0.5,0.5,0.5),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
  pie(c(nrow(func[[1]]),nrow(gene_length)-nrow(func[[1]])),
      labels=paste(c(paste0(dbs[1],":"),"not annotated:"),
                   c(nrow(func[[1]]),nrow(gene_length)-nrow(func[[1]]))),
      col =  c(brewer.pal(9,"Set1")[6],"grey"),clockwise = T,border=NA)
  dev.off()
}
