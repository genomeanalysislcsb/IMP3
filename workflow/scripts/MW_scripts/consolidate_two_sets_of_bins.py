#!/usr/bin/env python3
import sys
import os
from tabulate import tabulate

'''
This script takes in two folders containing bins from two binning methods, and their complementary 
contamination and completion scores derived from CheckM, then matches corresponding bins in the two
sets based on a minimum of 80% overlap (by length) of the bins, and finally, decides which of the two
bin versions is best in that particular bin.

Then it makes a new folder into which it puts the best version of each bin (changing the naming in
the process), and also makes a new .stats file which is consistent with the new bin folder.

Usage:
./script bin_folder_1 bin_folder_2 stats_file_1 stats_file_2 output_folder min_completion max_contamination
'''

# TODO: Mess, replace asap
def strip_tabulate_line(line):
    return line.strip().split()


def parse_tabulate_table_stats2list(file_path):
    with open(file_path) as f:
        data = []
        header = []
        for line in f:
            if line.startswith("-"):
                continue
            elif line.startswith("Bin"):
                if header:
                    print("Error: multiple headers in file")
                    raise ValueError
                header = strip_tabulate_line(line)
                data.append(header)
            else:
                data.append(strip_tabulate_line(line))
    return data


def load_bin_contigs(bin_files, bin_folder):
    """
    Loads the contig names and lengths for each bin file.

    Parameters:
    - bin_files: A list of filenames within the bin.
    - bin_folder: The folder path where the bin files are located.

    Returns:
    A dictionary where keys are bin filenames and values are dictionaries of contig names and lengths.
    """
    bins = {}
    for bin_file in bin_files:
        contig_len = 0
        contig_name = ""
        bins[bin_file] = {}
        with open(os.path.join(bin_folder, f"{bin_file}.fa")) as f:
            for line in f:
                if not line.startswith('>'):
                    contig_len += len(line.strip())
                else:
                    if contig_name != "":
                        bins[bin_file][contig_name] = contig_len
                        contig_len = 0
                    contig_name = line[1:-1]
        # Catch the last contig if the file doesn't end with a line staring on '>'
        if contig_name:
            bins[bin_file][contig_name] = contig_len
    return bins


# Get arguments
if len(sys.argv) != 8:
    print("Usage: ./script bin_folder_1 bin_folder_2 stats_file_1 stats_file_2 output_folder min_completion max_contamination")
    sys.exit(1)

bin_folder_1 = sys.argv[1]
bin_folder_2 = sys.argv[2]

stats_file_1 = sys.argv[3]
stats_file_2 = sys.argv[4]

output_folder = sys.argv[5]

min_completion = float(sys.argv[6])
max_contamination = float(sys.argv[7])


# load a list of good bins to save time (won't look at bad bins later on).
good_bins_1 = {}
good_bins_2 = {}

binner_1_bin_stats = parse_tabulate_table_stats2list(stats_file_1)

# Go through the stats file and set bad bins to None, skipping the header
for line in binner_1_bin_stats[1:]:
    bin_name = line[0]
    try:
        bin_comp = float(line[1])
        bin_cont = float(line[2])
    except ValueError:
        print(f"Error: {line}, file: {stats_file_1}, list: {binner_1_bin_stats}")
        raise ValueError
    if bin_comp > min_completion and bin_cont < max_contamination:
        good_bins_1[bin_name] = None

binner_2_bin_stats = parse_tabulate_table_stats2list(stats_file_2)
# Go through the stats file and set bad bins to None, skipping the header
for line in binner_2_bin_stats[1:]:
    bin_name = line[0]
    bin_comp = float(line[1])
    bin_cont = float(line[2])
    if bin_comp > min_completion and bin_cont < max_contamination:
        good_bins_2[bin_name] = None

# these are the dictionaries storing the contig names and lengths of each bin
bins_1 = load_bin_contigs(good_bins_1, bin_folder_1)
bins_2 = load_bin_contigs(good_bins_2, bin_folder_2)

print("Making all possible comparisons between the two bin sets, and recording total % identical length ...")
all_bin_pairs = {}
for bin_1 in good_bins_1:
    all_bin_pairs[bin_1] = {}
    for bin_2 in good_bins_2:
        # find identical contigs between bin_1 and bin_2
        match_1_length = 0
        match_2_length = 0
        mismatch_1_length = 0
        mismatch_2_length = 0

        for contig in bins_1[bin_1]:
            if contig in bins_2[bin_2]: match_1_length += bins_2[bin_2][contig]
            else: mismatch_1_length += bins_1[bin_1][contig]
        for contig in bins_2[bin_2]:
            if contig in bins_1[bin_1]: match_2_length += bins_1[bin_1][contig]
            else: mismatch_2_length += bins_2[bin_2][contig]

        # chose the highest % ID, depending on which bin is a subset of the other
        ratio_1 = 100 * match_1_length / (match_1_length + mismatch_1_length)
        ratio_2 = 100 * match_2_length / (match_2_length + mismatch_2_length)

        all_bin_pairs[bin_1][bin_2] = max([ratio_1, ratio_2])

print("Loading in completion and contamination scores of all the bins ...")
bins_1_stats = {}
bins_2_stats = {}
bins_1_summary = {}
bins_2_summary = {}

for i, bin_stats in enumerate(binner_1_bin_stats):
    if i == 0:
        bins_1_summary["header"] = bin_stats
        continue
    bin_name = bin_stats[0]
    if bin_name not in good_bins_1:
        continue
    bins_1_stats[bin_name] = (float(bin_stats[1]), float(bin_stats[2]))
    bins_1_summary[bin_name] = bin_stats[1:]

for i, bin_stats in enumerate(binner_2_bin_stats):
    if i == 0:
        bins_2_summary["header"] = bin_stats
        continue
    bin_name = bin_stats[0]
    if bin_name not in good_bins_2:
        continue
    bins_2_stats[bin_name] = (float(bin_stats[1]), float(bin_stats[2]))
    bins_2_summary[bin_name] = bin_stats[1:]

# go through all good bins and chose best ones
print("Going through the first group, pulling out identical bins from second group, and choosing the best version ...")
os.makedirs(output_folder, exist_ok=True)
selected_bins_stat_list = [bins_1_summary["header"]]
bins_2_matches = {}
bin_ct = 0
for bin_1 in all_bin_pairs:
    new_bin_name = f"bin.{bin_ct}"
    score = bins_1_stats[bin_1][0] - bins_1_stats[bin_1][1]*2
    found_better = False
    for bin_2 in all_bin_pairs[bin_1]:
        # check for sufficient overlap (80% bin length)
        if all_bin_pairs[bin_1][bin_2] < 80: continue
        bins_2_matches[bin_2] = None
        # check if this bin is better than original
        if (bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1]*2) > score:
            print(f"Found a better bin {bin_2} for {bin_1}.")
            selected_bins_stat_list.append([new_bin_name] + bins_2_summary[bin_2])
            cmd = f"cp {bin_folder_2}/{bin_2}.fa {output_folder}/{new_bin_name}.fa"
            os.system(cmd)
            found_better = True
    if not found_better:
        print(f"Keeping original bin {bin_1}.")
        selected_bins_stat_list.append([new_bin_name] + bins_1_summary[bin_1])
        cmd = f"cp {bin_folder_1}/{bin_1}.fa {output_folder}/{new_bin_name}.fa"
        os.system(cmd)
    bin_ct += 1   

print("Retrieving bins from second group that were not found in first group ...")
for bin_2 in bins_2_stats:
    new_bin_name = f"bin.{bin_ct}"
    if bins_2_stats[bin_2][0] < min_completion or bins_2_stats[bin_2][1] > max_contamination:
        continue
    if bin_2 in bins_2_matches:
        continue

    selected_bins_stat_list.append([new_bin_name] + bins_2_summary[bin_2])
    cmd = f"cp {bin_folder_2}/{bin_2}.fa {output_folder}/{new_bin_name}.fa"
    os.system(cmd)
    bin_ct += 1

# # Print all vars fo debugging
# print(f"good_bins_1: {good_bins_1}")
# print(f"good_bins_2: {good_bins_2}")
# print(f"bins_1: {bins_1}")
# print(f"bins_2: {bins_2}")
# print(f"all_bin_pairs: {all_bin_pairs}")
# print(f"bins_1_stats: {bins_1_stats}")
# print(f"bins_2_stats: {bins_2_stats}")
# print(f"bins_1_summary: {bins_1_summary}")
# print(f"bins_2_summary: {bins_2_summary}")
# print(f"selected_bins_stat_list: {selected_bins_stat_list}")
# print(f"bin_ct: {bin_ct}")
# print(f"bins_2_matches: {bins_2_matches}")
# print(f"output_folder: {output_folder}")

# raise Exception("Debugging")

# Write new stats file in tabulate format for consistency with the other scripts
# Sort the body based on the second column, in descending order
header = selected_bins_stat_list[0]
body = selected_bins_stat_list[1:]

sorted_body = sorted(body, key=lambda x: float(x[1]), reverse=True)
selected_bins_stat_list_sorted = [body] + sorted_body

out_table = tabulate(selected_bins_stat_list, headers="firstrow", tablefmt="simple")

print(f"Writing new stats file to {output_folder}/bin_stats.txt")

with open(f"{output_folder}.stats", "w") as f:
    f.write(out_table)

print("There were " + str(bin_ct) + " bins cherry-picked from the original sets!")
