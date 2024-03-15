#!/usr/bin/env python3
import sys
from tabulate import tabulate

# This script summarizes the statistics of each bin by parsing
col_list = ["Bin", "Completeness", "Contamination", "Contig_N50", "Genome_Size", "GC_Content", "binner"]

data = [col_list]

# Get input
if len(sys.argv) < 3 or len(sys.argv) > 4:
    print("Usage: ./summarize_checkm.py checkm_stats out_file [source_binner_stats_file]")
    sys.exit(1)
elif len(sys.argv) == 3:
    checkm_stats = sys.argv[1]
    out_file = sys.argv[2]
    binner_source_stats = None
    binner_source = checkm_stats.split(".")[0]
elif len(sys.argv) == 4:
    checkm_stats = sys.argv[1]
    out_file = sys.argv[2]
    binner_source_stats = sys.argv[3]

# If binner_source_stats is provided, get the binner source for each bin
if binner_source_stats:
    bin_source = {}
    # Get binning refiner combo/binning source name
    with open(binner_source_stats) as f:
        header = next(f)
        for line in f:
            try:
                cut = line.strip().split("\t")
                bin_source[cut[0]] = cut[-1]
            except IndexError:
                try:
                    cut = line.strip().split()
                    bin_source[cut[0]] = cut[-1]
                except IndexError:
                    print(f"Error splitting line with tabs or spaces: {line}")
                    raise IndexError

# Get CheckM stats
with open(checkm_stats) as file:
    header = next(file)
    for line in file:
        try:
            line_data = line.strip().split("\t")
        except IndexError:
            try:
                line_data = line.strip().split()
            except IndexError:
                print(f"Error splitting line with tabs or spaces: {line}")
                raise IndexError

        # Stats we need are: bin, comp, cont, n50, gen_size, gc
        bin_stats = [line_data[i] for i in [0, 1, 2, 6, 8, 9]]

        if binner_source_stats:
            bin_stats.append(bin_source.get(bin_stats[0], ''))   
        else:
            bin_stats.append(binner_source)
        
        data.append(bin_stats)

# Sort the data by completeness, ignoring the header but it it back later
# Separate the header and the data
header = data[0]
body = data[1:]

# Sort the body based on the second column, in descending order
sorted_body = sorted(body, key=lambda x: float(x[1]), reverse=True)

# Reassemble the data with the header at the beginning
sorted_data = [header] + sorted_body

# Now tabulate and print
print(tabulate(sorted_data, headers="firstrow", tablefmt="simple"))

# Write to file
with open(out_file, "w") as f:
    f.write(tabulate(sorted_data, headers="firstrow", tablefmt="simple"))
