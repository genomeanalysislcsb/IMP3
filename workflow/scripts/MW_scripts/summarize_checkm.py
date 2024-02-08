#!/usr/bin/env python3
import sys
import ast  # Import ast module to safely evaluate string literals as Python expressions

# This script summarizes the statistics of each bin by parsing 
# the checkm_folder/storage/bin_stats_ext.tsv file of the CheckM output

if len(sys.argv) == 3: 
    binner = sys.argv[2]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
elif len(sys.argv) == 4:
    source = {}
    with open(sys.argv[3]) as file:  # Use with statement for better file handling
        for line in file:
            cut = line.strip().split("\t")
            source[cut[0]] = cut[7]
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner")
else:
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize")

with open(sys.argv[1]) as file:  # Use with statement for better file handling
    for line in file:
        # Use ast.literal_eval for safer evaluation of string to dictionary
        dic = ast.literal_eval(line.strip().split("\t")[1])

        if "__" in dic["marker lineage"]: dic["marker lineage"] = dic["marker lineage"].split("__")[1]
        name = line.split("\t")[0]

        if len(sys.argv) == 3:    
            print("\t".join([name, str(dic["Completeness"])[:5],\
                str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
                dic["marker lineage"], str(dic["N50 (contigs)"]),\
                str(dic["Genome size"]), binner]))

        elif len(sys.argv) == 4:
            print("\t".join([name, str(dic["Completeness"])[:5],\
                str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
                dic["marker lineage"], str(dic["N50 (contigs)"]),\
                str(dic["Genome size"]), source.get(name, '')]))  # Use .get() to avoid KeyError

        else:
            print("\t".join([name, str(dic["Completeness"])[:5],\
                str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
                dic["marker lineage"], str(dic["N50 (contigs)"]),\
                str(dic["Genome size"])]))
