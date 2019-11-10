#!/bin/R
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# this script is used to merge mOTU and the taxonomic annotation of markergenes (of one COG) found by a BLAST search

#uses input from snakemake
#                  - COG, the COG ID of the genes of interest
#	           - besthitsfile, the best hit output of BLASTN "<COG>.bestPresentHits.tsv"  
#                  - name for the output
#also uses:
#			 - the mapping files "mOTU.v1.padded.motu.linkage.map", "mOTU.v1.map.txt" and "mOTU-LG.v1.annotations.txt" of the mOTU database
#output: a table with the query gene and the taxonomy of the best hit
# 

#written by Anna Heintz-Buschart (November 2019, based on the scripts from March 2019 and MuSt)

COG <- snakemake@wildcards[["COG"]]
blastoutput <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
ref_phylo <- snakemake@params[["refphy"]]
meta_phylo <- snakemake@params[["metaphy"]]

refphy <- read.delim(ref_phylo,stringsAsFactors=F)
metaphy <- read.delim(meta_phylo,stringsAsFactors=F)
phys <- data.frame("motuID"=c(refphy$ref.mOTU_v2_ID,metaphy$meta.mOTU_v2_ID),
                   rbind(refphy[,-c(1:2)],metaphy[,-1]),
                   stringsAsFactors=F)

bh <- read.delim(blastoutput,header=F,comment.char="#",stringsAsFactors=F)
colnames(bh) <- c("queryID", "subjectID", "percIdentity", "alignmentLength", "mismatches", "gapOpens", "q.start", "q.end", "s.start", "s.end", "evalue", "bitScore")

bh$subjectID <- sapply(bh$subjectID,function(x) unlist(strsplit(x,split=".",fixed=T))[1])
bh <- merge(bh,phys,by.x="subjectID",by.y="motuID",all.x=T)
write.table(bh,outfile,quote=F,sep="\t",row.names=F)


