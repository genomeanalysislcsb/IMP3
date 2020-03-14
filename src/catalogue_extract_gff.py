#!/usr/bin/env python3

import os
import sys

cluster_in = snakemake.input[0]
gff_in = snakemake.input[1] 
output = snakemake.output[0]

# to make gff from catalogue
# get rep genes from cluster_in
reps = {}
cluster_file = open(cluster_in, "r")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    rep = line.rstrip()
    reps[rep] = 1
cluster_file.close()
    
out_file = open(output, "w")
# extract lines from all samples' gffs for the rep genes
gff_file = open(gff_in, "r")
while 1:
    line = gff_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tabs = line.split("\t") #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
    attr = tabs[8]
    atts = attr.split(";")
    prot = [f.replace("ID=","") for f in atts if f.startswith("ID=")][0]
    if prot in reps:
        if tabs[2] == "CDS":
            parts = [f.replace("partial=","") for f in atts if f.startswith("partial=")][0]
            if tabs[6] == "-":
                partsr = parts[::-1]
                attr = attr.replace("partial=" + parts,"partial=" + partsr, 1)
        else:
            parts = "11" if "partial" in tabs[8] else "00"
            attr += ";partial="+parts
    # replace start and end with 1 and (end+1-start),
    # concatenate contig with gene name + sense + original coordinates for 1st column
    # make all sense +: 
        nc = tabs[0] + "_" + tabs[3] + "_" + tabs[4] + "_" + tabs[6]
        outl = [prot,tabs[1],"gene","1",str(int(tabs[4])+1-int(tabs[3])) ,tabs[5],"+",tabs[7],attr+";ori_cont="+nc]
        out_file.write("\t".join(outl) + "\n")
gff_file.close()
out_file.close()
