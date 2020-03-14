# this script is used to find the names of the genes (of one COG) in the mOTU padded database (Sunagawa et al 2013 Nat Methods) which represent the mOTUs present in a sample

#inputs: as arguments in the call:
#            - COG, the COG ID of the genes of interest
#			- presentMotus, the name of the file with the names of the present mOTUs
#			- output, the name of the output file
#also uses: 
#           - the mapping files "mOTU.v1.padded.motu.linkage.map" and "mOTU.v1.map.txt" of the mOTU database
#			- files with the names of the mOTU genes, extracted by COG from the mOTU database
#            - a prepared file with the names of the present mOTUs "mOTUs.present.all"  
#output: a list of genenames
# 

#written by Anna Heintz-Buschart (May 2108, based on the MuSt script from February 2015)

args<-commandArgs(TRUE)
COG <- args[1]
presentMotus <- args[2]
outfile <- args[3]

implib <- Sys.getenv("imp_lib")

phyl <- read.delim(paste0(implib,"/motus_linking/mOTU.v1.padded.motu.linkage.map"),stringsAsFactors=F)
map <- read.delim(paste0(implib,"/motus_linking/mOTU.v1.map.txt"),stringsAsFactors=F)
coord <- read.delim(paste0(implib,"/motus_linking/",COG,"/",COG,".mOTU.v1.padded.coord"),header=F,stringsAsFactors=F)
coord2 <- read.delim(paste0(implib,"/motus_linking/",COG,"/",COG,".mOTU.v1.padded.coord"),header=F,stringsAsFactors=F,sep=" ")
pres <- read.delim(presentMotus,header=F,stringsAsFactors=F)

phyl2 <- phyl[ phyl$mOTU.species.annotation %in% pres$V1,]
rm(phyl)
mp <- merge(map,phyl2,by.x="mOTU",by.y=1)
rm(list=c("map","phyl2"))
co <- gsub(" ",".",coord$V1)
rm(coord)
coord2$name <- co
mp2 <- merge(mp[mp$MarkerGene ==COG,],coord2[,c(1,4)],by.x="SequenceID",by.y=1)

write.table(mp2$name,outfile,col.names=F,quote=F,sep="\t",row.names=F)

