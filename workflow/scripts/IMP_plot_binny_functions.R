#!/bin/R

###################################################################################################
## This file contains all the functions required for all the plots
## 
###################################################################################################

## calculate coverage using the "traditional" method
get_coverage=function(reads_mapped, length, read_len){
    M <- reads_mapped
    L <- length
    R <- read_len

    # coverage*contig length/read length
    C <- (M*L)/R
    return(C)
}

####################################################################
## Contig based "RPKM" values, similar to used in
## Muller et al. (2014, Nat. Comm.), ## only here it
## is applied to
## every contig
contig_rpkm=function(reads_mapped, length){
    N <- sum(reads_mapped, na.rm=T)
    R <- reads_mapped
    L <- length

    # reads mapped/([length of contig]/1000)/([total reads]/10^6)
    R/(L/1000)/(N/10^6)
}

####################################################################
## Calculate variant density:
## i)  Traditional variats per kilo base
## ii) Normalized by contig rpkm (Muller et al., 2014)
var_density=function(variants, length, reads_mapped){
    V <- variants
    L <- length
    rpkmC <- contig_rpkm(reads_mapped, L)

    D <- (V/L)/1000/rpkmC
    return(D)
}

####################################################################
## Calculate gene density
gene_density=function(total_genes, length){
    G <- total_genes
    L <- length

    C <- (G)/(L/1000)
    return(C)
}
####################################################################
## Calculate coding density
coding_density=function(total_len_genes, length){
    G <- total_len_genes
    L <- length

    C <- G/L
    return(C)
}

####################################################################
## Calculate N50
get_n50=function(lengths){
    x=lengths
    x[cumsum(x) > sum(x, na.rm=T)/2][1]
}


####################################################################
## Get assembly statistics
get_stats=function(dat){
    contigs <- nrow(dat)
    N50 <- get_n50(dat$length)
    max_len <- max(dat$length)
    mean_len <- mean(dat$length)
    med_len <- median(dat$length)
    #MG_mapped <- sum(dat$MG_mapped)
    #MT_mapped <- sum(dat$MT_mapped)
    total_length <- sum(dat$length, na.rm=T)
    return(c(contigs,N50,max_len,mean_len,med_len,total_length))
}


####################################################################
## Function to scale data between 0 and 1 (for plotting)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

####################################################################
## Filter out outliers and set them as floor/ceiling value (min/max)
outliers=function(z, dist){
    z[is.infinite(z)] <- max(z[is.finite(z)])
    z[which(z > mean(z) + dist*sd(z))] = round(mean(z) + dist*sd(z))
    z[which(z < mean(z) - dist*sd(z))] = round(mean(z) - dist*sd(z))
    return(z)
}

####################################################################
## Function to name files
name_plot=function(name){
    filename <- paste(out_dir, name, sep='/')
    return(filename)
}

####################################################################
# Function to obtain number of character (char) occurences
# within string
countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}

####################################################################
# Function to output filename without the path
get_file_name <- function(file_path){
    n <- countCharOccurrences('/', as.character(file_path))
    filename <- str_split_fixed(file_path, "/", n+1)[n+1]
    return(filename)
}

####################################################################
# Function to obtain filtering information (contained in file names)
filtering <- function(file){
    n <- countCharOccurrences("\\.", as.character(file))
    filter <- str_split_fixed(file, "\\.", n+1)[n]
    return(filter)
}

## Function for white background theme with no axes
theme_nothing <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
	    panel.border = element_rect(fill="NA", colour="NA"),
            legend.text = element_text(size=base_size*1.5, hjust=1.2, vjust=1), 
            legend.title = element_text(size=base_size*2,face="bold",hjust=0),
	    legend.key = element_rect(fill="NA", colour="white"),
	    legend.key.width = unit(1.2, "cm"),
	    legend.key.height = unit(1, "cm"),
	    legend.position = "right",
	    legend.direction = "vertical",
	    legend.box = "vertical",
            line = element_blank(),
            rect             = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}

## Function for black background theme with no axes
theme_black <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
	    panel.background = element_rect(fill="black", colour="black"),
	    panel.border = element_rect(fill="NA", colour="NA"),
            legend.text = element_text(size=base_size*1.5), 
            legend.title = element_text(size=base_size*2,face="bold",hjust=0),
	    legend.key = element_rect(fill="NA", colour="white"),
	    legend.key.height = unit(1, "cm"),
	    legend.key.width = unit(1.2, "cm"),
	    legend.position = "right",
	    legend.direction = "vertical",
	    legend.box = "vertical",
            line = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}

## Function for gray background theme with no axes
theme_gray <- function(base_size = 12, base_family = "Helvetica")
  {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
      theme(
	    panel.background = element_rect(fill="gray25", colour="gray25"),
            line             = element_blank(),
            axis.ticks.margin = unit(0, "lines"),
	    axis.text.x=element_blank(),
	    axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
           )
}

## Function for beanplot
mgmt.beanplot <- function(MG.dat, MT.dat, ylabel = expression(log[10]*~"count")) {
	var1 <-log10(c(MG.dat, MT.dat))
	var1[is.infinite(var1)]=NA
	var2 <- c(rep("MG", length(MG.dat)), rep("MT", length(MT.dat)))
	var <- data.frame(var1, var2)
	
	par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
	beanplot(var1 ~ var2, 
		 data= var,  
		 side = "both",
		 log = "auto",
	what=c(1,1,1,0), 
	border = NA, 
	col = list("blue", c("red", "white")),
	bw = "nrd0", 
	ylab = ylabel)
	legend("bottomleft", fill = c("blue", "red"), legend = c("MG", "MT"))
}
## Set maximum value for plots based on standard deviation
set_max_sd=function(x, dist){
    max_val <- mean(x, na.rm=T) + dist*sd(x, na.rm=T)
    return(max_val)
}

## Set maximum value for plots based on percentage
set_max_perc=function(x, percentage){
    max_val <- max(x[is.finite(x)], na.rm=T)*(percentage/100)
}

## Special string for length
log10len <- expression(bold(atop("Contig length", paste("(",log[10], bp, ")"))))

## Metagenomic and metatranscriptomic labels
mgmt_labs <- c("metagenomic","metatranscriptomic")

## Set plot height and width

plotWidth = 1200
plotHeight = 1000

## binny functions
find.cutoff <- function(data,k=2,proba=0.5) {
  model <- try(normalmixEM(x=data, k=k,maxrestarts=50))
  if(class(model)!="try-error"){
    i <- which.min(model$mu)
    f <- function(x) {
      proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
    }
    lower <- min(model$mu)
    upper <- max(model$mu)
    if(f(lower)*f(upper)>0){
      return(median(data))
    }else{
      return(uniroot(f=f, lower=lower, upper=upper)$root)  # Careful with division by zero if changing lower and upper
    }
  }else{
    return(median(data))
  }
}

muClus <- function(clusterName,cRes,resFile,recDepth){
  covcol <- which(colnames(contigInfo)=="MG_depth")
  dInfo <- contigInfo[cRes$cluster==clusterName,]
  if(dip.test(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))$p.value<0.05){
    cuto <- find.cutoff(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))
    write.table(t(c(clusterName,cuto)),resFile,append=T,row.names=F,col.names=F,quote=F,sep="\t")
    subs1 <- dInfo[log10(1+dInfo[,covcol])<cuto,]
    num1 <- length(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni1 <- length(unique(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni1==0) {
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else if(num1/uni1<=1.2){
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else {
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(cRes$cluster[cRes$contig %in% subs1$contig],1,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs1$contig]),cRes,resFile,recDepth+1)
      print(recDepth)
    }
    subs2 <- dInfo[log10(1+dInfo[,covcol])>=cuto,]
    num2 <- length(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni2 <- length(unique(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni2==0) {
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else if(num2/uni2<=1.2){
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else {
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(cRes$cluster[cRes$contig %in% subs2$contig],2,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs2$contig]),cRes,resFile,recDepth+1)
      print(recDepth)
    }
  } else {
    cRes$cluster[cRes$contig %in% dInfo$contig] <- gsub("^.","B",cRes$cluster[cRes$contig %in% dInfo$contig])
  }
  return(cRes)
}


