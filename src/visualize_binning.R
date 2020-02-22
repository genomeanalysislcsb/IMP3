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
library(alluvial)
library(data.table)

###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input[["stats"]] #"Stats/all_stats.Rdata"
bin_memberships     <- snakemake@input[["bins"]]  # Binning/*/scaffold2bin.tsv
mag_memberships     <- snakemake@input[["dastool"]]  # Binning/selected_DASTool_scaffolds2bin.txt
outplot             <- snakemake@output # "Visualization/Assembly_length_hist.png"
#                                         Visualization/Assembly_cov_hist.png"
                                        #"Visualization/Assembly_gc_hist.png"
                                         #"Visualization/Assembly_cov_gc_scatter.png"

###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

binners <- unique(binTab$binner)

###################################################################################################
## plot
###################################################################################################
print("Plotting alluvial diagram.")
dir.create("Visualization", showWarnings = FALSE)

plotfile <- outplot[[1]]
cont2bin <- list()
for(binner in binners){
  cont2bin[[binner]] <- read.delim(grep(binner,bin_memberships,value=T),stringsAsFactors = F,header=F)
  colnames(cont2bin[[binner]]) <- c("contig","bin")
}
cont2bin[["DASTool"]] <- read.delim(mag_memberships,stringsAsFactors = F,header=F)
colnames(cont2bin[["DASTool"]]) <- c("contig","bin")

cont2binM <- rbindlist(lapply(1:length(cont2bin),function(x) data.frame(cont2bin[[x]],
                                                                       "tool"=names(cont2bin)[x],
                                                                       stringsAsFactors = F)))
tools <- levels(as.factor(cont2binM$tool))
cont2binM <- apply(tapply(cont2bin$bin,list(cont2bin$contig,cont2bin$tool),paste),1,function(x)paste(x,sep="|",collapse="|"))
cont2binM <- merge(cont2binM,contig_length,by.x=0,by.y=1,all.x=T)
cont2binM <- aggregate(cont2binM$length,list(cont2binM$x),sum)
cont2binM$Group.1 <- as.character(cont2binM$Group.1)
cont2binM <- data.frame(t(sapply(cont2binM$Group.1,function(x) unlist(strsplit(x,split="|",fixed=T)))),
                       "freq"=cont2binM$x,stringsAsFactors = F)
colnames(cont2binM)[1:length(tools)] <- tools
#prettier names
if("MaxBin" %in% colnames(cont2binM)){
  cont2binM$MaxBin <- gsub("maxbin_res.","",gsub(".fasta","",cont2binM$MaxBin))
  cont2binM$DASTool <- gsub("maxbin_res.","",gsub(".fasta","",cont2binM$DASTool))
}
all_col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*sapply(cont2binM$DASTool,function(x) {
  if(x %in% binTab$selected_by_DASTool) binTab$binScore.DASTool[binTab$selected_by_DASTool==x] else 0
})]
png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
par("mar"=c(0.2,0.2,0.2,0.2),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
alluvial(cont2binM[,c(binners,"DASTool")],freq = cont2binM$freq,
         col=all_col,
         border=NA)
dev.off()


if("mg" %in% types){
  png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
  par("mar"=c(13,4.5,0.5,0.5),tcl="-0.3",mgp=c(3.2,0.5,0),lwd=0.5)
  cont2binC <- merge(cont2bin[["DASTool"]],depthAssStat[["mg"]],by=1,all.x=T)
  cont2binC[is.na(cont2binC)] <- 0
  plotnames <- sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                      function(x) paste0(x,"\n",gsub(".+;(.__[^;]+).*","\\1",binTab$classification[binTab$selected_by_DASTool==x])))
  box_col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                                                                           function(x) binTab$binScore.DASTool[binTab$selected_by_DASTool==x] )]
  zeroC <- sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                  function(x) length(which(cont2binC$bin[cont2binC$bin==x]==0)))
    a<-boxplot(cont2binC$depth[cont2binC$depth>0]~cont2binC$bin[cont2binC$depth>0],las=2,
          axes=F,log="y",ylab="",pch=16,cex=0.5,border=box_col)
  text(1:length(unique(cont2binC$bin[cont2binC$depth>0])),
       rep(min(cont2binC$depth[cont2binC$depth>0])*0.9,length(unique(cont2binC$bin[cont2binC$depth>0]))),
       labels = paste(zeroC,"0s"))
  mtext("contig coverage depths",2,3.2,cex=9/7)
  axis(2,las=1,lwd=0.5)
  axis(1,at=1:length(unique(cont2binC$bin[cont2binC$depth>0])),
       labels=plotnames,las=2,lwd=0.5)
  box()
  dev.off()
}

if("mt" %in% types){
  png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
  par("mar"=c(13,4.5,0.5,0.5),tcl="-0.3",mgp=c(3.2,0.5,0),lwd=0.5)
  tmpMT <- merge(c2g[c2g$contig %in% cont2bin[["DASTool"]]$contig,],
                 depthGeneStat[["mt"]],by="gene",all.x=T)
  tmpMT[is.na(tmpMT)] <- 0
  cont2binT <- aggregate(tmpMT$depth,list(tmpMT$cont),function(x) length(which(x>0))/length(x))
  colnames(cont2binT) <- c("contig","percExpr")
  cont2binT <- merge(cont2bin[["DASTool"]],cont2binT,by="contig")
  totE <- cont2binT <- merge(tmpMT,cont2binT,by="contig")
  totE <- aggregate(totE$depth,list(totE$bin),function(x) length(which(x>0))/length(x))
  plotnames <- sapply(levels(as.factor(cont2binT$bin)),
                      function(x) paste0(x,"\n",gsub(".+;(.__[^;]+).*","\\1",binTab$classification[binTab$selected_by_DASTool==x])))
  box_col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*sapply(levels(as.factor(cont2binT$bin)),
                                                                                  function(x) binTab$binScore.DASTool[binTab$selected_by_DASTool==x] )]
  a<-boxplot(100*cont2binT$percExpr~cont2binT$bin,las=2,
             axes=F,ylab="",pch=16,cex=0.5,border=box_col,
             ylim=c(0,100),outline=F)
  text(1:length(unique(cont2binT$bin)),
       rep(5,length(unique(cont2binT$bin))),
       labels = paste(round(100*totE$x),"%"))
  axis(2,at=5,labels="overall\nexpressed",las=2,tick = F,lwd=0.4)
  mtext("% expressed genes per contig",2,3.2,cex=9/7)
  axis(2,las=1,lwd=0.5)
  axis(1,at=1:length(unique(cont2binT$bin)),
       labels=plotnames,las=2,lwd=0.5)
  box()
  dev.off()
}

if("binny" %in% binners){
  png(plotfile,width=17,height=17,units = "cm",res=300,pointsize=7)
  par(mar=c(0.5,0.5,0.5,0.5))
  "Binning/binny/mgmt.vizbin.with-contig-names.points"
  plotlim <- 1.05*max(abs(c(min(contigInf[,coordCol[1]]),min(contigInf[,coordCol[2]]),
                            max(contigInf[,coordCol[1]]),max(contigInf[,coordCol[2]]))))
spec <- colorRampPalette(brewer.pal(9,"Spectral")[9:1])
if(class(contigInf[,colCol])=="character"){
  plotcol <- spec(length(unique(contigInf[contigInf[,colCol]!="unknown",colCol])))[as.numeric(as.factor(contigInf[contigInf[,colCol]!="unknown",colCol]))]
}else{
  plotcol <- spec(length(unique(contigInf[!is.na(contigInf[,colCol]),colCol])))[as.numeric(as.factor(contigInf[!is.na(contigInf[,colCol]),colCol]))]
}
plot(contigInf[contigInf[,colCol]!="unknown"&!is.na(contigInf[,colCol]),coordCol],pch=16,cex=0.3,las=1,col=plotcol,
     ylim=c(-plotlim,plotlim),xlim=c(-plotlim,plotlim),ann=F,axes=F)
points(contigInf[contigInf[,colCol]=="unknown"|is.na(contigInf[,colCol]),coordCol],pch=16,cex=0.3,las=1,col="grey80")
box(lwd=2)
if(essHelper){
  essHelp <- unique(unlist(sapply(c("GrpE","TIGR00115","TIGR01030"),function(x)grep(x,contigInf$essentialGene))))
  if(length(essHelp)>0){
    ehCh <- sapply(contigInf$essentialGene[essHelp],function(x) if(any(grepl("GrpE",x))) 8 else if(any(grepl("TIGR00115",x))) 4 else 3)
    points(contigInf[essHelp,coordCol],pch=ehCh,cex=0.8,col="grey35")
  }
}
if(class(contigInf[,colCol])=="character"){
  legend("topleft",legend=paste(c(unique(contigInf[contigInf[,colCol]!="unknown",colCol]),"unknown"),sapply(c(unique(contigInf[contigInf[,colCol]!="unknown",colCol]),"unknown"),function(x)(length(which(contigInf[,colCol]==x))))),
         col=c(unique(plotcol),"grey80"),pch=16,cex=0.4,bty="n",y.intersp=0.8,title=paste(nrow(contigInf),"contigs"))
}else{
  legend("topleft",legend=paste(colnames(contigInf)[colCol],sort(unique(contigInf[!is.na(contigInf[,colCol]),colCol]))[seq(from=1,to=length(unique(contigInf[!is.na(contigInf[,colCol]),colCol])),length.out=10)]),
         col=spec(length(unique(contigInf[!is.na(contigInf[,colCol]),colCol])))[seq(from=1,to=length(unique(contigInf[!is.na(contigInf[,colCol]),colCol])),length.out=10)],pch=16,cex=0.4,bty="n",y.intersp=0.8,title=paste(nrow(contigInf),"contigs"))
}
if(essHelper){
  legend("topright",legend=paste(c("GrpE","TIGR00115","TIGR01030"),"-",
                                 sapply(c("GrpE","TIGR00115","TIGR01030"),function(x)(length(which(contigInf$essentialGene==x)))),"times"),
         col="grey35",pch=c(8,4,3),cex=0.4,bty="n",y.intersp=0.8,title="essential Gene markers")
}
mtext(sampleID,3,-1.5,font=2,cex=0.8)
dev.copy(pdf,paste(gsub("-","",sub("20","",Sys.Date())),"_",name,".pdf",sep=""),
         pointsize=9,width=7,height=7)
dev.off()


}
# vizbin scatter plot with coverage
#circos prep