#!/usr/bin/env python3

import os
import sys

cluster_in = snakemake.input[0]
gff_in = snakemake.input[1:] 
output = snakemake.output

# to make gff from catalogue
# get rep genes from cluster_in
reps = {}
cluster_file = open(cluster_in, "r")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    line = line.rstrip()
    rep = line.split("\t")[0]
    reps[rep] = 1
cluster_file.close()
    
out_file = open(output, "w")
# extract lines from all samples' gffs for the rep genes
for gff in gff_in:
    gff_file = open(gff, "r")
    while 1:
        line = gff_file.readline()
        if line == "":
            break
        line = line.rstrip()
        tabs = line.split("\t") #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
        atts = tabs[8].split(";")
        prot = [f.replace("ID=","") for f in atts if f.startswith("ID=")][0]
        if prot in reps:
    # replace start and end with 1 and (end+1-start),
    # concatenate contig with gene name + sense + original coordinates for 1st column
    # make all sense +: 
            nc = tabs[0] + "_" + prot + "_" + tabs[3] + "_" + tabs[4] + "_" + tabs[6]
            outl = [nc,tabs[1],tabs[2],1,tabs[4]+1-tabs[3] ,tabs[5],"+",tabs[7],tabs[8]]
            out_file.write("\t".join(outl) + "\n")
    gff_file.close()
out_file.close()
