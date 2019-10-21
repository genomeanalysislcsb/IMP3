# this script is used to merge mOTU and the taxonomic annotation of markergenes (of one COG) found by a BLAST search

#inputs: as argument in the call 
#                  - COG, the COG ID of the genes of interest
#				   - besthitsfile, the best hit output of BLASTN "<COG>.bestPresentHits.tsv"  
#				  - name for the oupt
#also uses:
#			 - the mapping files "mOTU.v1.padded.motu.linkage.map", "mOTU.v1.map.txt" and "mOTU-LG.v1.annotations.txt" of the mOTU database
#output: a table with the query gene and the taxonomy of the best hit
# 

#written by Anna Heintz-Buschart (May 2018, based on the MuSt script from Maarch 2015)

args<-commandArgs(TRUE)
COG <- args[1]
blastoutput <- args[2]
outfile <- args[3]

implib <- Sys.getenv("imp_lib")

phyl <- read.delim(paste0(implib,"/motus_linking/mOTU.v1.padded.motu.linkage.map"),stringsAsFactors=F)
map <- read.delim(paste0(implib,"/motus_linking/mOTU.v1.map.txt"),stringsAsFactors=F)
mapPhyl <- read.delim(paste0(implib,"/motus_linking/mOTU-LG.v1.annotations.txt"),stringsAsFactors=F)

bh <- read.delim(blastoutput,header=F,comment.char="#",stringsAsFactors=F)
colnames(bh) <- c("queryID", "subjectID", "percIdentity", "alignmentLength", "mismatches", "gapOpens", "q.start", "q.end", "s.start", "s.end", "evalue", "bitScore")

bh$subjectID <- sapply(bh$subjectID,function(x) paste(unlist(strsplit(x,split=".",fixed=T))[-c(length(unlist(strsplit(x,split=".",fixed=T)))-1:0)],sep=".",collapse="."))
bh <- merge(bh,map,by.x=2,by.y=1,all.x=T)
rm(map)
bh <- merge(bh,phyl[,c(1,11)],by.x="mOTU",by.y=1,all.x=T)
rm(phyl)
bh <- merge(bh,mapPhyl,by.x=ncol(bh),by.y=1,all.x=T)
rm(mapPhyl)
write.table(bh,outfile,quote=F,sep="\t",row.names=F)


