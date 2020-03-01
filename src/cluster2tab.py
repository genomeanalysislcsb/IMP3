#!/usr/bin/env python3

import os
import sys
import pandas as pd

cluster_in = snakemake.input
output = snakemake.output
samples_in = snakemake.params.sample
samples = samples_in[samples_in.analysis == 1]

# get rep genes from cluster_in
# count number of members per cluster + number of members per cluster per sample
cluster_file = open(cluster_in, "r")
out_file = open(output, "w")
curr_clus = ""
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    line = line.rstrip()
    rep = line.split("\t")[0]
    mem = line.split("\t")[1]
    sam = samples.loc[samples.prokka_prefix == mem.split("_")[0],"sample"]
    if rep != curr_clus:
        out_file.write("\t".join([str(rep[s]) for s in samples.sample.append("all")]) + "\n")
        curr_clus = rep
        reps = dict(zip(samples.sample.append("all"),[0]*(1+len(samples.sample))))
    reps[rep]["all"] += 1
    reps[rep][sam] += 1
cluster_file.close()
out_file.write("\t".join([str(rep[s]) for s in samples.sample.append("all")]) + "\n")
out_file.close()
    
