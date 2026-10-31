#!/usr/bin/env python3
import sys, os

# Usage: ./script.py bins.stats binsFolder outFolder

binsStats = sys.argv[1]
binsFolder = sys.argv[2]
outFolder = sys.argv[3]
remove_dups = False
if len(sys.argv) > 4:
    if sys.argv[4] == "remove":
        remove_dups = True

# load in bin completion and contamination scores
print("Loading in bin completion and contamination scores...")
bin_scores = {}
with open(binsStats) as file:
    header = next(file)
    seperator_line = next(file)
    for line in file:
        bin_data = line.strip().split()
        bin_name = bin_data[0]
        bin_comp = float(bin_data[1])
        bin_cont = float(bin_data[2])
        bin_N50 = float(bin_data[3])
        # Score = completion - 2 * contamination + 0.0000000001 * contig_N50
        try:
            score = bin_comp - 2 * bin_cont + 0.0000000001 * bin_N50
        except IndexError:
            print(f"Error: {line}, file: {binsStats}")
            raise IndexError
        bin_scores[bin_name] = score

# load in contigs in each bin
print("Loading in contigs in each bin...")
contig_mapping = {}
for bin_file in os.listdir(binsFolder):
    bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
    with open(os.path.join(binsFolder, bin_file)) as file:
        for line in file:
            if line[0] != ">": continue
            contig = line[1:-1]
            if contig not in contig_mapping: contig_mapping[contig] = bin_name
            else:
                if remove_dups:
                    contig_mapping[contig] = None
                elif bin_scores[bin_name] > bin_scores[contig_mapping[contig]]:
                    contig_mapping[contig] = bin_name

# go over the bin files again and make a new dereplicated version of each bin file
print("Making a new dereplicated version of each bin file")
os.makedirs(outFolder, exist_ok=True)
for bin_file in os.listdir(binsFolder):
    bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
    with open(os.path.join(outFolder, bin_file), 'w') as out:
        at_least_one = False
        store = False
        for line in open(os.path.join(binsFolder, bin_file)):
            if line[0] == ">":
                contig = line[1:-1]
                if contig_mapping.get(contig) == bin_name:
                    at_least_one = True
                    store = True
                else:
                    store = False
            if store: out.write(line)
        if not at_least_one:
            os.remove(os.path.join(outFolder, bin_file))
