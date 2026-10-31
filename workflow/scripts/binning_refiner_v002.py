#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

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

# Version 1.4.3, modified by Oskar Hickl, 2024-02-01

import os
import shutil
import argparse
import gzip
from Bio import SeqIO
from datetime import datetime
import sys


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)

    folder_list_sorte = sorted(folder_list)
    return folder_list_sorte


def is_gzipped(file_path):
    """Check if a file is gzipped by looking at its magic bytes."""
    try:
        with open(file_path, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except:
        return False


def open_fasta(filename, mode='r'):
    """Open a FASTA file, automatically handling gzipped files."""
    if is_gzipped(filename):
        return gzip.open(filename, f"{mode}t")  # Open in text mode
    else:
        return open(filename, mode)


def plot_sankey(data_file, output_file, width=800, height=600):
    """Create a sankey plot using pysankey2."""
    try:
        # Try importing pysankey2
        try:
            from pysankey2 import Sankey
            import pandas as pd
            import matplotlib.pyplot as plt
            
            # Read data
            df = pd.read_csv(data_file)
            
            # If the CSV has only one column with comma-separated values
            if len(df.columns) == 1:
                header = df.columns[0].split(',')
                # Parse CSV properly
                df = pd.read_csv(data_file, skiprows=1, header=None)
                # Split the values by comma if they're in a single column
                if len(df.columns) == 1:
                    df = df[0].str.split(',', expand=True)
                
                # Rename the columns
                if len(header) >= 3:
                    layer_cols = header[:2]  # Take first two elements for layer names
                    df.columns = layer_cols + [f'layer{i+3}' for i in range(len(df.columns)-2)]
                    
                    # Convert weight column to numeric
                    df[df.columns[2]] = pd.to_numeric(df[df.columns[2]], errors='coerce')
                    
                    # Create Sankey diagram
                    sankey = Sankey(df, layers=df.columns[:2], values=df.columns[2], colorMode="global")
                    fig, ax = sankey.plot(figsize=(width/100, height/100))
                    
                    # Save figure
                    output_png = output_file.replace('.html', '.png')
                    fig.savefig(output_png, dpi=100, bbox_inches='tight')
                    plt.close(fig)
                    
                    print(f"Sankey plot created at {output_png}")
                    return True
            else:
                # Handle properly formatted CSVs directly
                cols = df.columns.tolist()
                if len(cols) >= 3:
                    sankey = Sankey(df, layers=cols[:2], values=cols[2], colorMode="global")
                    fig, ax = sankey.plot(figsize=(width/100, height/100))
                    
                    # Save figure
                    output_png = output_file.replace('.html', '.png')
                    fig.savefig(output_png, dpi=100, bbox_inches='tight')
                    plt.close(fig)
                    
                    print(f"Sankey plot created at {output_png}")
                    return True
            
            print("CSV format not as expected, skipping sankey plot")
            return False
                
        except ImportError:
            print("pysankey2 module not found. Install with: pip install git+https://github.com/SZJShuffle/pySankey2.git")
            return False
            
    except Exception as e:
        print(f"Error creating sankey plot: {str(e)}")
        return False


######################################## CONFIGURATION #######################################

parser = argparse.ArgumentParser()

parser.add_argument('-i',       required=True,  help='input bin folder')
parser.add_argument('-o',       required=False, default='Bin_refiner',      help='output folder, default: Bin_refiner')
parser.add_argument('-p',       required=False, default='Refined',      help='output prefix, default: Refined')
parser.add_argument('-m',       required=False, default=512,  type=int, help='minimal size (Kbp) of refined bin, default: 512')
parser.add_argument('-plot',    required=False, default=False,          action="store_true", help='visualize refinement with Sankey plot')
parser.add_argument('-x',       required=False, default=False,          help='the width of sankey plot')
parser.add_argument('-y',       required=False, default=False,          help='the height of sankey plot')
parser.add_argument('-q',       required=False, default=False,          help='silent progress report')
parser.add_argument('-gz',      required=False, default=True,          action="store_true", help='compress output with gzip')

args = vars(parser.parse_args())

input_bin_folder =      args['i']
output_dir =            args['o']
output_bin_prefix =     args['p']
minBin_size_Kbp =       args['m']
plot_sankey_flag =      args['plot']
sankey_plot_width =     args['x']
sankey_plot_height =    args['y']
keep_quiet =            args['q']
compress_output =       args['gz']


if input_bin_folder[-1] == '/':
    input_bin_folder = input_bin_folder[:-1]


################################### define output filename ###################################

pwd_report_file_sources_and_length =  '%s/%s_sources_and_length.txt'    % (output_dir, output_bin_prefix)
pwd_report_file_contigs =             '%s/%s_contigs.txt'               % (output_dir, output_bin_prefix)
pwd_report_file_for_sankey =          '%s/%s_sankey.csv'                % (output_dir, output_bin_prefix)
pwd_plot_file_for_sankey =            '%s/%s_sankey.html'               % (output_dir, output_bin_prefix)
pwd_output_bin_folder =               '%s/%s'                           % (output_dir, 'Refined')

# Create output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(pwd_output_bin_folder, exist_ok=True)

time_format = '[%Y-%m-%d %H:%M:%S]'


################################### precheck of input files ##################################

# get input_bin_subfolders
input_bin_subfolders = get_no_hidden_folder_list(input_bin_folder)

for bin_subfolder in input_bin_subfolders:
    # get bin file list in each input bin subfolder
    pwd_bin_subfolder = '%s/%s' % (input_bin_folder, bin_subfolder)
    bin_file_list = get_no_hidden_folder_list(pwd_bin_subfolder)

    bin_ext_list = set()
    for each_bin in bin_file_list:
        each_bin_ext = each_bin.split('.')[-1]
        bin_ext_list.add(each_bin_ext)

    if len(bin_ext_list) > 1:
        print('Program exited, please make sure all bins within %s folder have the same extension' % bin_subfolder)
        exit()


##################################### refine input bins ######################################

# progress report
if keep_quiet is False:
    print('%s Detected %s input bin folders' % ((datetime.now().strftime(time_format)), len(input_bin_subfolders)))

ctg_length_dict = {}
ctg_to_bin_dict = {}
for bin_subfolder in input_bin_subfolders:

    # get bin file list in each input bin subfolder
    pwd_bin_subfolder = '%s/%s' % (input_bin_folder, bin_subfolder)
    bin_file_list = get_no_hidden_folder_list(pwd_bin_subfolder)

    # progress report
    if keep_quiet is False:
        print('%s Read in %s %s bins' % ((datetime.now().strftime(time_format)), len(bin_file_list), bin_subfolder))

    for each_bin in bin_file_list:
        pwd_each_bin = '%s/%s/%s' % (input_bin_folder, bin_subfolder, each_bin)

        try:
            with open_fasta(pwd_each_bin) as f:
                for each_seq in SeqIO.parse(f, 'fasta'):
                    each_seq_id = str(each_seq.id)
                    each_seq_len = len(each_seq.seq)

                    # store contig length into dict
                    if each_seq_id not in ctg_length_dict:
                        ctg_length_dict[each_seq_id] = each_seq_len

                    # store contig to bin info into dict
                    if each_seq_id not in ctg_to_bin_dict:
                        ctg_to_bin_dict[each_seq_id] = [each_bin]
                    else:
                        ctg_to_bin_dict[each_seq_id].append(each_bin)
        except Exception as e:
            print(f"Error reading {pwd_each_bin}: {str(e)}")
            continue

# progress report
if keep_quiet is False:
    print('%s Refine input bins' % (datetime.now().strftime(time_format)))


# only keep contigs existed in all bin sets
ctg_to_bin_dict_shared = {}
for each_ctg in ctg_to_bin_dict:
    if len(ctg_to_bin_dict[each_ctg]) == len(input_bin_subfolders):
        ctg_to_bin_dict_shared[each_ctg] = '___'.join(ctg_to_bin_dict[each_ctg])


# get concatenated_bins_to_ctg_dict
concatenated_bins_to_ctg_dict = {}
concatenated_bins_length_dict = {}
for each_shared_ctg in ctg_to_bin_dict_shared:

    concatenated_bins = ctg_to_bin_dict_shared[each_shared_ctg]

    if concatenated_bins not in concatenated_bins_to_ctg_dict:
        concatenated_bins_to_ctg_dict[concatenated_bins] = [each_shared_ctg]
        concatenated_bins_length_dict[concatenated_bins] = ctg_length_dict[each_shared_ctg]
    else:
        concatenated_bins_to_ctg_dict[concatenated_bins].append(each_shared_ctg)
        concatenated_bins_length_dict[concatenated_bins] += ctg_length_dict[each_shared_ctg]


# remove short length
concatenated_bins_to_ctg_dict_short_removed = {}
concatenated_bins_length_dict_short_removed = {}
maximum_length = 0
for concatenated_bins in concatenated_bins_length_dict:
    concatenated_bins_length = concatenated_bins_length_dict[concatenated_bins]
    if concatenated_bins_length >= minBin_size_Kbp * 1024:
        concatenated_bins_to_ctg_dict_short_removed[concatenated_bins] = concatenated_bins_to_ctg_dict[concatenated_bins]
        concatenated_bins_length_dict_short_removed[concatenated_bins] = concatenated_bins_length_dict[concatenated_bins]

        if concatenated_bins_length_dict[concatenated_bins] > maximum_length:
            maximum_length = concatenated_bins_length_dict[concatenated_bins]


# add leading zero to length
concatenated_bins_length_dict_short_removed_str = {}
for each_length in concatenated_bins_length_dict_short_removed:
    concatenated_bins_length_dict_short_removed_str[each_length] = ("{:0%sd}" % len(str(maximum_length))).format(concatenated_bins_length_dict_short_removed[each_length])


concatenated_bins_list = []
for each_concatenated_bin in concatenated_bins_length_dict_short_removed_str:
    concatenated_bins_list.append('%s___%s' % (concatenated_bins_length_dict_short_removed_str[each_concatenated_bin], each_concatenated_bin))


# sort concatenated_bins
concatenated_bins_list_sorted = sorted(concatenated_bins_list, reverse=True)


# progress report
if keep_quiet is False:
    print('%s Got %s refined bins with size larger than %s Kbp' % ((datetime.now().strftime(time_format)), len(concatenated_bins_list_sorted), minBin_size_Kbp))


# progress report
if keep_quiet is False:
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Extract sequences for refined bins'))


# extract sequences and write out report file
pwd_report_file_sources_and_length_handle = open(pwd_report_file_sources_and_length, 'w')
pwd_report_file_sources_and_length_handle.write('Refined_bin\tSize(Kbp)\tSource\n')

pwd_report_file_contigs_handle = open(pwd_report_file_contigs, 'w')
pwd_report_file_contigs_handle.write('Refined_bin\tContigs\n')

pwd_report_file_for_sankey_handle = open(pwd_report_file_for_sankey, 'w')
pwd_report_file_for_sankey_handle.write('C1,C2,Length_Kbp\n')

refined_bin_index = 1
for refined_bin in concatenated_bins_list_sorted:

    refined_bin_split = refined_bin.split('___')
    refined_bin_length = concatenated_bins_length_dict_short_removed['___'.join(refined_bin_split[1:])]
    refined_bin_length_Kbp = float("{0:.2f}".format(refined_bin_length/1024))
    refined_bin_source = refined_bin_split[1:]
    refined_bin_ctgs = concatenated_bins_to_ctg_dict['___'.join(refined_bin_split[1:])]
    refined_bin_id = '%s_%s' % (output_bin_prefix, refined_bin_index)
    
    # Set output path based on compression choice
    if compress_output:
        pwd_refined_bin = '%s/%s_%s.fasta.gz' % (pwd_output_bin_folder, output_bin_prefix, refined_bin_index)
    else:
        pwd_refined_bin = '%s/%s_%s.fasta' % (pwd_output_bin_folder, output_bin_prefix, refined_bin_index)
    
    pwd_ctg_source = '%s/%s/%s' % (input_bin_folder, input_bin_subfolders[0], refined_bin_split[1])

    # write out sources_and_length file
    pwd_report_file_sources_and_length_handle.write('%s\t%s\t%s\n' % (refined_bin_id, refined_bin_length_Kbp, ','.join(refined_bin_source)))

    # write out contig file
    pwd_report_file_contigs_handle.write('%s\t%s\n' % (refined_bin_id, ','.join(sorted(refined_bin_ctgs))))

    # write out file for sankey plot
    n = 0
    while n <= len(refined_bin_source) - 2:
        pwd_report_file_for_sankey_handle.write('%s\n' % ','.join([refined_bin_source[n], refined_bin_source[n+1], str(refined_bin_length_Kbp)]))
        n += 1

    # extract sequences for refined bin
    try:
        # Open output file with appropriate compression
        if compress_output:
            pwd_refined_bin_handle = gzip.open(pwd_refined_bin, 'wt')
        else:
            pwd_refined_bin_handle = open(pwd_refined_bin, 'w')
            
        with open_fasta(pwd_ctg_source) as f:
            for each_source_ctg in SeqIO.parse(f, 'fasta'):
                if str(each_source_ctg.id) in refined_bin_ctgs:
                    SeqIO.write(each_source_ctg, pwd_refined_bin_handle, 'fasta')
        
        pwd_refined_bin_handle.close()
    except Exception as e:
        print(f"Error extracting sequences from {pwd_ctg_source}: {str(e)}")

    refined_bin_index += 1

pwd_report_file_sources_and_length_handle.close()
pwd_report_file_contigs_handle.close()
pwd_report_file_for_sankey_handle.close()


# progress report
if plot_sankey_flag:
    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), 'Creating Sankey plot'))
    
    width = sankey_plot_width if sankey_plot_width else len(input_bin_subfolders) * 200
    height = sankey_plot_height if sankey_plot_height else len(concatenated_bins_list_sorted) * 12
    
    # Try using pysankey
    plot_success = plot_sankey(pwd_report_file_for_sankey, pwd_plot_file_for_sankey, width, height)
    
    if not plot_success and not keep_quiet:
        print('%s %s' % ((datetime.now().strftime(time_format)), 'Sankey plot creation failed or skipped'))


# progress report
if keep_quiet is False:
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Done!'))