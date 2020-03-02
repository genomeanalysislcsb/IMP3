#!/usr/bin/env python3

import os
import sys
import pandas as pd

cluster_in = snakemake.input[0]
output = snakemake.output[0]
samples_in = snakemake.params.sample
samples = samples_in[samples_in.genes == 1]


# get rep genes from cluster_in
# count number of members per cluster + number of members per cluster per sample
cluster_file = open(cluster_in, "r")
out_file = open(output, "w")
curr_clus = ""
all_s = samples["sample"].tolist()
all_s.append("all")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    line = line.rstrip()
    rep = line.split("\t")[0]
    if curr_clus== "":
        curr_clus = rep
        reps = dict(zip(all_s,[0]*len(all_s)))
        out_file.write("gene\t" + "\t".join(all_s) + "\n")
    mem = line.split("\t")[1]
    sam = samples.loc[samples.prokka_prefix == mem.split("_")[0],"sample"].iat[0]
    if rep != curr_clus:
        out_file.write(curr_clus + "\t" + "\t".join([str(reps[s]) for s in all_s]) + "\n")
        curr_clus = rep
        reps = dict(zip(all_s,[0]*len(all_s)))
    reps["all"] += 1
    reps[sam] += 1
cluster_file.close()
out_file.write(curr_clus + "\t" + "\t".join([str(reps[s]) for s in all_s]) + "\n")
out_file.close()
    
