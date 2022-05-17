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

makeCircosFiles <- function(contigs,genes,col,circos_files){
 dir.create(paste0(circos_files,"/data"), showWarnings = FALSE) 
  kt <- data.frame("chr"=rep("chr",nrow(contigs)),
                   "ps"=rep("-",nrow(contigs)),
                   "ID"=contigs$contig,"LABEL"=contigs$contig,
                   "START"=rep(0,nrow(contigs)),"END"=contigs$length-1,
                   "COLOR"=rep(col,nrow(contigs)),stringsAsFactors=F)
  kt <- kt[order(kt$END,decreasing=T),]
  write.table(kt,paste0(circos_files,"/data/contigs.ideogram",sep=""),row.names=F,col.names=F,quote=F)
  kt <- data.frame("chr"=genes$contig[genes$strand=="+"&genes$feature=="CDS"],
                   "START"=genes$start[genes$strand=="+"&genes$feature=="CDS"],
                   "END"=genes$end[genes$strand=="+"&genes$feature=="CDS"],
                   "COLOR"=paste("color=",ifelse(genes$completeness[genes$strand=="+"&genes$feature=="CDS"]=="complete",
                                                 "complete","incomplete"),sep=""),
                   stringsAsFactors=F)
  write.table(kt,paste0(circos_files,"/data/features.forward.txt",sep=""),row.names=F,col.names=F,quote=F)
  kt <- data.frame("chr"=genes$contig[genes$strand=="-"&genes$feature=="CDS"],
                   "START"=genes$start[genes$strand=="-"&genes$feature=="CDS"],
                   "END"=genes$end[genes$strand=="-"&genes$feature=="CDS"],
                   "COLOR"=paste("color=",ifelse(genes$completeness[genes$strand=="-"&genes$feature=="CDS"]=="complete",
                                                 "complete","incomplete"),sep=""),
                   stringsAsFactors=F)
  write.table(kt,paste0(circos_files,"/data/features.reverse.txt",sep=""),row.names=F,col.names=F,quote=F)
  if(length(genes$contig[genes$strand=="+"&genes$feature=="rRNA"])>0){
    kt <- data.frame("chr"=genes$contig[genes$strand=="+"&genes$feature=="rRNA"],
                     "START"=genes$start[genes$strand=="+"&genes$feature=="rRNA"],
                     "END"=genes$end[genes$strand=="+"&genes$feature=="rRNA"],
                     "COLOR"=paste("color=",ifelse(genes$completeness[genes$strand=="+"&genes$feature=="rRNA"]=="complete",
                                                   "completerrna","incompleterrna"),sep=""),
                     stringsAsFactors=F)
  } else {
    kt <- ""
  }
  write.table(kt,paste0(circos_files,"/data/rRNAs.forward.txt",sep=""),row.names=F,col.names=F,quote=F)
  if(length(genes$contig[genes$strand=="-"&genes$feature=="rRNA"])>0){
    kt <- data.frame("chr"=genes$contig[genes$strand=="-"&genes$feature=="rRNA"],
                     "START"=genes$start[genes$strand=="-"&genes$feature=="rRNA"],
                     "END"=genes$end[genes$strand=="-"&genes$feature=="rRNA"],
                     "COLOR"=paste("color=",ifelse(genes$completeness[genes$strand=="-"&genes$feature=="rRNA"]=="complete",
                                                   "completerrna","incompleterrna"),sep=""),
                     stringsAsFactors=F)
  } else {
    kt <- ""
  }
  write.table(kt,paste0(circos_files,"/data/rRNAs.reverse.txt",sep=""),row.names=F,col.names=F,quote=F)
  if("mg" %in% colnames(genes)){
    kt <- data.frame("chr"=genes$contig,"START"=genes$start,"END"=genes$end,
                     "value"=log10(0.001+genes$mg),stringsAsFactors=F)
  }else {
    kt <- ""
  }
  write.table(kt,paste0(circos_files,"/data/features.DNAcovGenes.txt",sep=""),row.names=F,col.names=F,quote=F)
  
  if("mt" %in% colnames(genes)){
    if(any(genes$mt>0)){
      rnaval <- log10(0.001+genes$mt)
    }else{
      rnaval <- rep(0,nrow(genes))
    }
  kt <- data.frame("chr"=genes$contig,"START"=genes$start,"END"=genes$end,
                   "value"=rnaval,stringsAsFactors=F)
  }else {
    kt <- ""
  }
  write.table(kt,paste0(circos_files,"/data/features.RNAcovGenes.txt",sep=""),row.names=F,col.names=F,quote=F)
  if("essential" %in% colnames(genes)){
    kt <- data.frame("chr"=genes$contig[genes$essential!=""],
                     "START"=genes$start[genes$essential!=""],
                     "END"=genes$end[genes$essential!=""])
    tabEss <- table(genes$essentialGene[genes$essential!=""])
    lt <- data.frame("chr1","START1","END1","chr2","START2","END2",stringsAsFactors=F)[F,]
    if(any(tabEss>1)){
      dupEss <- names(tabEss[tabEss>1])
      for(de in dupEss){
        dupGe <- genes[genes$essential==de,c("contig","start","end")]
        for(i in 1:(nrow(dupGe)-1)){
          for(j in 2:nrow(dupGe)){
            lt[nrow(lt)+1,] <- c(dupGe[i,1],
                                 dupGe[i,2:3],
                                 dupGe[j,1],
                                 dupGe[j,2:3])
          }
        }
      }
    }
  }else {
    kt <- ""
    lt <- ""
  }
  write.table(kt,paste0(circos_files,"/data/essentialGenes.highlights.txt",sep=""),row.names=F,col.names=F,quote=F)
  write.table(lt,paste0(circos_files,"/data/essentialGenes.links.txt",sep=""),row.names=F,col.names=F,quote=F)
}

gffGetID <- function(gff){
 gff$gene <- gsub("ID=([^;]+);.+","\\1",gff$attribute)
  return(gff)
}

gffAddComplete <- function(gff){
  gff$completeness <- sapply(1:nrow(gff), function(x){
    if(gff$feature[x]=="CDS"){
      if(gsub(".+partial=([01]{2});.*","\\1",paste0(gff$attribute[x],";"))=="00") "complete" else "incomplete"
    }else{
      if(!grepl("partial",gff$attribute[x])) "complete" else "incomplete"
    }
  })
  return(gff)                        
}

gffAddFunc <- function(gff,func){
  gff[[func]] <- sapply(1:nrow(gff), function(x){
    if(gff$feature[x]=="CDS" & grepl(paste0(".+;",func,"=([^;]+)"),gff$attribute[x])){
      gsub(paste0(".+;",func,"=([^;]+);.*"),"\\1",paste0(gff$attribute[x],";"))
    }else{
      ""
    }
  })
  return(gff) 
}

min0 <- function(x){
  min(x[x>0])
}

spec <- colorRampPalette(brewer.pal(9,"Spectral")[9:1])
ygb <- colorRampPalette(brewer.pal(8,"YlGnBu")[2:8])

###################################################################################################
## Files from snakemake
###################################################################################################

stats_image         <- snakemake@input[["stats"]] #"Stats/all_stats.Rdata"
bin_memberships     <- snakemake@input[["bins"]]  # Binning/*/scaffold2bin.tsv
mag_memberships     <- snakemake@input[["dastool"]]  # Binning/selected_DASTool_scaffolds2bin.txt
if(file.exists("Binning/binny"))  vizbin_file         <- snakemake@input[["vizbin"]] #"Binning/binny/ASS.vizbin.with-contig-names.points",
outplot             <- snakemake@output # "Visualization/Binning_alluvial.png"
                                        #"Visualization/Binning_{type}_depth.png"
                                       #"Visualization/Binning_vizbin.png" 
circos_files       <- "Visualization/"

###################################################################################################
## Load image 
###################################################################################################
print("Loading data.")
load(stats_image)

if(!"binTab" %in% ls()){
 for(o in unlist(outplot)){
      system(paste0("touch ",o))
    }
    print("No binning stats. Aborting visualization.")
}else{

binners <- unique(binTab$binner)
names(outplot)[1:(1+length(types))] <- c("alluvial",types)
if(file.exists("Binning/binny")) names(outplot)[length(outplot)] <- "vizbin"

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
cont2binM <- apply(tapply(cont2binM$bin,list(cont2binM$contig,cont2binM$tool),paste),1,function(x)paste(x,sep="|",collapse="|"))
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
  if(x %in% binTab$selected_by_DASTool) binTab$binScore.DASTool[which(binTab$selected_by_DASTool==x)] else 0
})]
png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
par("mar"=c(0.2,0.2,0.2,0.2),tcl="-0.3",mgp=c(1.7,0.5,0),lwd=0.5)
alluvial(cont2binM[,c(binners,"DASTool")],freq = cont2binM$freq,
         col=all_col,
         border=NA)
dev.off()

if("depthGeneStat" %in% ls()){
for(bin in unique(cont2bin[["DASTool"]]$bin)){
  circ_genes <- gffGetID(gene_gff[gene_gff$contig %in% cont2bin[["DASTool"]]$contig[cont2bin[["DASTool"]]$bin==bin],])
  circ_genes <- gffAddComplete(circ_genes)
  if("essential" %in% dbs | "binny" %in% binners){
    circ_genes <- gffAddFunc(circ_genes,"essential")
  }
  if("mg" %in% types){
    circ_genes <- merge(circ_genes,
                        depthGeneStat[["mg"]],by="gene",all.x=T)
    colnames(circ_genes)[ncol(circ_genes)] <- "mg"
    circ_genes$mg[is.na(circ_genes$mg)] <- 0
  }
  if("mt" %in% types){
    circ_genes <- merge(circ_genes,
                        depthGeneStat[["mt"]],by="gene",all.x=T)
    colnames(circ_genes)[ncol(circ_genes)] <- "mt"
    circ_genes$mt[is.na(circ_genes$mt)] <- 0
  }
  col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*binTab$binScore.DASTool[binTab$selected_by_DASTool==bin]]
  makeCircosFiles(contig_length[contig_length$contig %in% cont2bin[["DASTool"]]$contig[cont2bin[["DASTool"]]$bin==bin],],
                   circ_genes,col,circos_files)
}}
if("mg" %in% types){
  plotfile <- outplot[["mg"]]
  png(plotfile,width=17,height=12,units = "cm",res=300,pointsize=7)
  par("mar"=c(13,4.5,0.5,0.5),tcl="-0.3",mgp=c(3.2,0.5,0),lwd=0.5)
  cont2binC <- merge(cont2bin[["DASTool"]],depthAssStat[["mg"]],by=1,all.x=T)
  cont2binC[is.na(cont2binC)] <- 0
  plotnames <- sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                      function(x) paste0(x,"\n",gsub(".+;(.__[^;]+).*","\\1",binTab$classification[which(binTab$selected_by_DASTool==x)])))
  box_col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                                                                           function(x) binTab$binScore.DASTool[which(binTab$selected_by_DASTool==x)] )]
  zeroC <- sapply(levels(as.factor(cont2binC$bin[cont2binC$depth>0])),
                  function(x) length(which(cont2binC$bin[cont2binC$bin==x]==0)))
    a<-boxplot(cont2binC$depth[cont2binC$depth>0]~cont2binC$bin[cont2binC$depth>0],las=2,
          axes=F,log="y",ylab="",pch=16,cex=0.5,border=box_col)
  text(1:length(unique(cont2binC$bin[cont2binC$depth>0])),
       rep(min(cont2binC$depth[cont2binC$depth>0])*0.9,length(unique(cont2binC$bin[cont2binC$depth>0]))),
       labels = paste(zeroC,"0s"),cex=0.5)
  mtext("contig coverage depths",2,3.2,cex=9/7)
  axis(2,las=1,lwd=0.5)
  axis(1,at=1:length(unique(cont2binC$bin[cont2binC$depth>0])),
       labels=plotnames,las=2,lwd=0.5,cex=0.5)
  box()
  dev.off()
}
if("depthGeneStat" %in% ls()){
plotfile <- outplot[["mt"]]
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
                      function(x) paste0(x,"\n",gsub(".+;(.__[^;]+).*","\\1",binTab$classification[which(binTab$selected_by_DASTool==x)])))
  box_col <- colorRampPalette(brewer.pal(11,"Spectral"))(256)[256:1][1+255*sapply(levels(as.factor(cont2binT$bin)),
                                                                                  function(x) binTab$binScore.DASTool[which(binTab$selected_by_DASTool==x)] )]
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
}
if(file.exists("Binning/binny")){
plotfile <- outplot[["vizbin"]]
#if("binny" %in% binners){
  png(plotfile,width=17,height=17,units = "cm",res=300,pointsize=7)
  par(mar=c(0.5,0.5,0.5,0.5))
  vizbin <- read.delim(vizbin_file,stringsAsFactors = F,header=F)
  colnames(vizbin) <- c("contig","x","y")
  vizcont <- merge(vizbin,depthAssStat[["mg"]],by="contig",all.x=T)
  vizcont[is.na(vizcont)] <- 0
  plotlim <- 1.05*max(abs(c(min(vizcont$x),min(vizcont$y),
                            max(vizcont$x),max(vizcont$y))))
  plotcol <- ygb(length(unique(vizcont$depth)))[as.numeric(as.factor(vizcont$depth))]
  plot(vizcont[,c("x","y")],pch=16,cex=0.3,las=1,col=plotcol,
       ylim=c(-plotlim,plotlim),xlim=c(-plotlim,plotlim),ann=F,axes=F)
  box(lwd=0.5)
  legend("topleft",
         legend=paste("coverage depth",round(sort(unique(vizcont$depth))[seq(from=1,to=length(unique(vizcont$depth)),
                                                        length.out=10)])),
         col=ygb(length(unique(vizcont$depth)))[seq(from=1,to=length(unique(vizcont$depth)),
                                                 length.out=10)],
         pch=16,cex=0.8,bty="n",y.intersp=0.8,title=paste(nrow(vizcont),"contigs"))
  dev.off()
}
}
