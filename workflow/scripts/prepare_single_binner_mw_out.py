#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 8 11:58:32 2022

@author: oskar.hickl
"""

import sys
import shutil


def parse_mw_stats_file(mw_stats_file, mw_comp_thresh, mw_cont_thresh, mw_stats_header=True,
                        write_filtered_stats_file=True):
    tmp_stats_file = mw_stats_file + '.tmp'
    good_bin_list = []
    with open(mw_stats_file, 'r') as f:
        if mw_stats_header:
            next(f)
        for line in f:
            line_data = line.split('\t')
            bin = line_data[0]
            comp = float(line_data[1])
            cont = float(line_data[2])
            if comp > mw_comp_thresh and cont < mw_cont_thresh:
                good_bin_list.append(bin)
                if write_filtered_stats_file:
                    with open(tmp_stats_file, 'a') as f2:
                        f2.write(line)
    if write_filtered_stats_file:
        shutil.move(tmp_stats_file, mw_stats_file)
    return good_bin_list


def prepare_single_binner_mw_out(mw_stats_file, binner_bin_dir, comp_thresh, cont_thresh, mw_dir, mw_stats_header=True):
    good_bins = parse_mw_stats_file(mw_stats_file, comp_thresh, cont_thresh, mw_stats_header,
                                    write_filtered_stats_file=True)
    # Copy bins, write bins.contigs file
    with open(f'{mw_dir}/bins.contigs', 'w') as bc:
        for bin in good_bins:
            shutil.copy2(f'{binner_bin_dir}/{bin}.fasta', f'{mw_dir}/bins/{bin}.fasta')
            with open(f'{mw_dir}/bins/{bin}.fasta', 'r') as bf:
                for line in bf:
                    if line.startswith('>'):
                        line_contig = line.strip('>\n')
                        bc.write(f'{line_contig}\t{bin}\n')


mw_stats = sys.argv[1]
binner_bin_dir = sys.argv[2]
comp_thresh = float(sys.argv[3])
cont_thresh = float(sys.argv[4])
mw_dir = sys.argv[5]
header = eval(sys.argv[6])

prepare_single_binner_mw_out(mw_stats, binner_bin_dir, comp_thresh, cont_thresh, mw_dir, header)





