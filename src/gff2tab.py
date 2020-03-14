#!/usr/bin/env python3

import os
import sys

gff = snakemake.input[0]
output = snakemake.output[0]
dbs = snakemake.params.dbs

gff_file = open(gff, "r")    
out_file = open(output, "w")
# get length, functions and completeness
out_file.write("\t".join(["gene","length","completeness"]+dbs) + "\n")
while 1:
    line = gff_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tabs = line.split("\t") #0: contig 1:source 2:feature 3:start 4:end 5:. 6:sense 7:score 8:attribute
    atts = tabs[8].split(";")
    prot = [f.replace("ID=","") for f in atts if f.startswith("ID=")][0]
    part = [f.replace("partial=","") for f in atts if f.startswith("partial=")][0]
    if int(part) == 0:
        compl = "complete"
    else:
        compl = "incomplete"
    func=[]
    for db in dbs:
        func.append(";".join([f.replace(db + "=","") for f in atts if f.startswith(db + "=")]))
    out_file.write("\t".join([prot,tabs[4],compl] + func) + "\n")
gff_file.close()
out_file.close()
