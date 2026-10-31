#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com
# t.thomas@unsw.edu.au

# Binning_refiner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Binning_refiner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# metaWRAP author notes:
# I thank the original creator of this script! This is a great idea! To make
# this script more usable as part of the metaWRAP binning pipeline, I
# removed unnecessary visual aaspects of the original Bin_refiner script
# and made it python2 compatible.

# Check out the original program: https://github.com/songweizhi/Binning_refiner
# And the publication: https://www.ncbi.nlm.nih.gov/pubmed/28186226

# IMP note:
# Basic port to python3 was done by Oskar Hickl, 2024

import os
import glob
import shutil
import argparse
from time import sleep
from sys import stdout
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument('-1', required=True, help='first bin folder name')
parser.add_argument('-2', required=True, help='second bin folder name')
parser.add_argument('-3', required=False, help='third bin folder name')
parser.add_argument('-o', required=True, help='output folder name')
parser.add_argument('-ms', required=False, default=524288, type=int, help='(optional) minimum size for refined bins, default = 524288 (0.5Mbp)')

args = vars(parser.parse_args())
output_dir = args['o']
if output_dir[-1] == '/':
    output_dir = output_dir[:-1]

input_bin_folder_1 = args['1']
if input_bin_folder_1[-1] == '/':
    input_bin_folder_1 = input_bin_folder_1[:-1]

input_bin_folder_2 = args['2']
if input_bin_folder_2[-1] == '/':
    input_bin_folder_2 = input_bin_folder_2[:-1]

if args['3'] is not None:
    input_bin_folder_3 = args['3']
    if input_bin_folder_3[-1] == '/':
        input_bin_folder_3 = input_bin_folder_3[:-1]

bin_size_cutoff = args['ms']
bin_size_cutoff_MB = float("{0:.2f}".format(bin_size_cutoff / (1024 * 1024)))

input_bin_folder_list = []
if args['3'] is None:
    print('Specified 2 input bin sets: -1 {} -2 {}'.format(input_bin_folder_1, input_bin_folder_2))
    input_bin_folder_list = [input_bin_folder_1, input_bin_folder_2]
else:
    print('Specified 3 input bin sets: -1 {} -2 {} -3 {}'.format(input_bin_folder_1, input_bin_folder_2, input_bin_folder_3))
    input_bin_folder_list = [input_bin_folder_1, input_bin_folder_2, input_bin_folder_3]

wd = os.getcwd()
print('Working directory: {}'.format(wd))

output_folder = output_dir
pwd_output_folder = '{}/{}'.format(wd, output_folder)
os.mkdir(pwd_output_folder)
print('Output folder: {}'.format(pwd_output_folder))


folder_bins_dict = {}
all_input_bins_list = []
all_input_bins_number_list = []
for bin_folder in input_bin_folder_list:
    bins_files = '{}/{}/*.fa'.format(wd, bin_folder)
    bin_folder_bins = [os.path.basename(file_name) for file_name in glob.glob(bins_files)]
    all_input_bins_list.append(bin_folder_bins)
    all_input_bins_number_list.append(len(bin_folder_bins))
    folder_bins_dict[bin_folder] = bin_folder_bins
    if not bin_folder_bins:
        print('No input bin detected from {} folder, please double-check!'.format(bin_folder))
        exit()

    bin_folder_bins_ext_list = [os.path.splitext(bin)[1] for bin in bin_folder_bins]
    bin_folder_bins_ext_list_uniq = list(set(bin_folder_bins_ext_list))
    if len(bin_folder_bins_ext_list_uniq) > 1:
        print('Different file extensions were found from {} bins, please use the same extension (fa, fas, or fasta) for all bins in the same folder.'.format(bin_folder))
        exit()

if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)

combined_all_bins_file = '{}/{}/combined_all_bins.fasta'.format(wd, output_folder)
separator = '__'
for each_folder in input_bin_folder_list:
    sleep(1)
    print('Add folder/bin name to contig name for {} bins'.format(each_folder))
    os.mkdir('{}/{}/{}_new'.format(wd, output_folder, each_folder))
    each_folder_bins = folder_bins_dict[each_folder]
    for each_bin in each_folder_bins:
        bin_file_name, bin_file_ext = os.path.splitext(each_bin)
        each_bin_content = SeqIO.parse('{}/{}/{}'.format(wd, each_folder, each_bin), 'fasta')
        with open('{}/{}/{}_new/{}_{}.fasta'.format(wd, output_folder, each_folder, each_folder, bin_file_name), 'w') as new:
            for each_contig in each_bin_content:
                each_contig_new_id = '{}{}{}{}{}'.format(each_folder, separator, bin_file_name, separator, each_contig.id)
                each_contig.id = each_contig_new_id
                each_contig.description = ''
                SeqIO.write(each_contig, new, 'fasta')

    os.system('cat {}/{}/{}_new/*.fasta > {}/{}/combined_{}_bins.fa'.format(wd, output_folder, each_folder, wd, output_folder, each_folder))
    os.system('rm -r {}/{}/{}_new'.format(wd, output_folder, each_folder))

print('Combine all bins together')
if len(input_bin_folder_list) == 2:
    pwd_combined_folder_1_bins = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_1)
    pwd_combined_folder_2_bins = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_2)
    os.system('cat {} {} > {}'.format(pwd_combined_folder_1_bins, pwd_combined_folder_2_bins, combined_all_bins_file))
elif len(input_bin_folder_list) == 3:
    pwd_combined_folder_1_bins = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_1)
    pwd_combined_folder_2_bins = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_2)
    pwd_combined_folder_3_bins = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_3)
    os.system('cat {} {} {} > {}'.format(pwd_combined_folder_1_bins, pwd_combined_folder_2_bins, pwd_combined_folder_3_bins, combined_all_bins_file))

combined_all_bins = SeqIO.parse(combined_all_bins_file, 'fasta')
contig_bin_dict = {}
contig_length_dict = {}
for each in combined_all_bins:
    each_id_split = each.id.split(separator)
    folder_name = each_id_split[0]
    bin_name = each_id_split[1]
    contig_id = each_id_split[2]
    length = len(each.seq)
    if contig_id not in contig_bin_dict:
        contig_bin_dict[contig_id] = ['{}{}{}'.format(folder_name, separator, bin_name)]
        contig_length_dict[contig_id] = length
    else:
        contig_bin_dict[contig_id].append('{}{}{}'.format(folder_name, separator, bin_name))

contig_assignments_file = '{}/{}/contig_assignments.txt'.format(wd, output_folder)
with open(contig_assignments_file, 'w') as contig_assignments:
    for each in contig_bin_dict:
        if len(contig_bin_dict[each]) == len(input_bin_folder_list):
            contig_assignments.write('{}\t{}\t{}\n'.format('\t'.join(contig_bin_dict[each]), each, contig_length_dict[each]))

contig_assignments_file_sorted = '{}/{}/contig_assignments_sorted.txt'.format(wd, output_folder)
contig_assignments_file_sorted_one_line = '{}/{}/contig_assignments_sorted_one_line.txt'.format(wd, output_folder)
os.system('cat {} | sort > {}'.format(contig_assignments_file, contig_assignments_file_sorted))

with open(contig_assignments_file_sorted) as contig_assignments_sorted, open(contig_assignments_file_sorted_one_line, 'w') as contig_assignments_sorted_one_line:
    current_match = ''
    current_match_contigs = []
    current_length_total = 0
    n = 1
    for each in contig_assignments_sorted:
        each_split = each.strip().split('\t')
        current_contig = each_split[-2]
        current_length = int(each_split[-1])
        matched_bins = '\t'.join(each_split[:-2])
        if current_match == '':
            current_match = matched_bins
            current_match_contigs.append(current_contig)
            current_length_total += current_length
        elif current_match == matched_bins:
            current_match_contigs.append(current_contig)
            current_length_total += current_length
        elif current_match != matched_bins:
            refined_bin_name = 'refined_bin{}'.format(n)
            if current_length_total >= bin_size_cutoff:
                contig_assignments_sorted_one_line.write('Refined_{}\t{}\t{}bp\t{}\n'.format(n, current_match, current_length_total, '\t'.join(current_match_contigs)))
                n += 1
            current_match = matched_bins
            current_match_contigs = [current_contig]
            current_length_total = current_length

    if current_length_total >= bin_size_cutoff:
        contig_assignments_sorted_one_line.write('Refined_{}\t{}\t{}bp\t{}\n'.format(n, current_match, current_length_total, '\t'.join(current_match_contigs)))

print('The number of refined bins: {}'.format(n))

print('Exporting refined bins...')
separated_1 = '{}/{}/Refined_bins_sources_and_length.txt'.format(wd, output_folder)
separated_2 = '{}/{}/Refined_bins_contigs.txt'.format(wd, output_folder)
googlevis_input_file = '{}/{}/GoogleVis_Sankey_{}Mbp.csv'.format(wd, output_folder, bin_size_cutoff_MB)
os.mkdir('{}/{}/Refined'.format(wd, output_folder))
with open(contig_assignments_file_sorted_one_line) as refined_bins, open(googlevis_input_file, 'w') as googlevis_input_handle, open(separated_1, 'w') as separated_1_handle, open(separated_2, 'w') as separated_2_handle:
    googlevis_input_handle.write('C1,C2,Length (Mbp)\n')
    for each_refined_bin in refined_bins:
        each_refined_bin_split = each_refined_bin.strip().split('\t')
        each_refined_bin_name = each_refined_bin_split[0]
        each_refined_bin_length = 0
        each_refined_bin_contig = []
        if len(input_bin_folder_list) == 2:
            each_refined_bin_source = each_refined_bin_split[1:3]
            each_refined_bin_length = int(each_refined_bin_split[3][:-2])
            each_refined_bin_contig = each_refined_bin_split[4:]
            separated_1_handle.write('{}\t{}bp\t{}\n'.format(each_refined_bin_name, each_refined_bin_length, '\t'.join(each_refined_bin_source)))
            separated_2_handle.write('{}\n{}\n'.format(each_refined_bin_name, '\t'.join(each_refined_bin_contig)))

        elif len(input_bin_folder_list) == 3:
            each_refined_bin_source = each_refined_bin_split[1:4]
            each_refined_bin_length = int(each_refined_bin_split[4][:-2])
            each_refined_bin_contig = each_refined_bin_split[5:]
            separated_1_handle.write('{}\t{}bp\t{}\n'.format(each_refined_bin_name, each_refined_bin_length, '\t'.join(each_refined_bin_source)))
            separated_2_handle.write('{}\n{}\n'.format(each_refined_bin_name, '\t'.join(each_refined_bin_contig)))
        each_refined_bin_length_mbp = float("{0:.2f}".format(each_refined_bin_length / (1024 * 1024)))
        m = 0
        while m < len(each_refined_bin_source)-1:
            googlevis_input_handle.write('{},{},{}\n'.format(each_refined_bin_source[m], each_refined_bin_source[m+1], each_refined_bin_length_mbp))
            m += 1

        stdout.write('Extracting refined bin: {}.fasta\n'.format(each_refined_bin_name))
        refined_bin_file = '{}/{}/Refined/{}.fasta'.format(wd, output_folder, each_refined_bin_name)
        with open(refined_bin_file, 'w') as refined_bin_handle:
            input_contigs_file = '{}/{}/combined_{}_bins.fa'.format(wd, output_folder, input_bin_folder_1)
            input_contigs = SeqIO.parse(input_contigs_file, 'fasta')
            for each_input_contig in input_contigs:
                each_input_contig_id = each_input_contig.id.split(separator)[-1]
                if each_input_contig_id in each_refined_bin_contig:
                    each_input_contig.id = each_input_contig_id
                    each_input_contig.description = ''
                    SeqIO.write(each_input_contig, refined_bin_handle, 'fasta')

# Check if all output dirs and files were created
out_files = [googlevis_input_file, separated_1, separated_2]
for out_file in out_files:
    if not os.path.isfile(out_file):
        print('Error: {} was not created!'.format(out_file))
        exit()
    else:
        print('Created: {}'.format(out_file))

# Remove temporary files
print('\nDeleting temporary files')
# os.system('rm {}'.format(contig_assignments_file))
# os.system('rm {}'.format(combined_all_bins_file))
# os.system('rm {}'.format(contig_assignments_file_sorted))
# os.system('rm {}'.format(contig_assignments_file_sorted_one_line))
print('\nAll done!')
