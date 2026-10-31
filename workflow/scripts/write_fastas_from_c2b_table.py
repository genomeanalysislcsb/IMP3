#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 11:34:52 2022

@author: oskar.hickl
"""

import sys


def load_fasta(FASTA,DICT):
    """

    Parameters
    ----------
    FASTA : STR
        Path to input fasta.
    DICT : Dict
       Dict to put fasta headers as keys and sequences as list into.

    Returns
    -------
    Dict with headers as keys and sequences as list.

    """
    with open(FASTA, 'r') as f:
        current_header = None
        header_line = False
        for line in f:
            if line.startswith(">"):
                header_line = True
            else:
                header_line = False
            line = line.strip(">\n\t ")
            # Get header
            if header_line:
                # Transform sublists of sequences from previous header to one list.
                if current_header:
                    DICT[current_header] = "".join(DICT[current_header])
                DICT[line] = []
                current_header = line
            # Get sequences
            else:
                DICT[current_header].append(line)
        # Transform sublists of sequences from last header to one list.
        DICT[current_header] = "".join(DICT[current_header])
    return DICT


assembly = sys.argv[1]
binner_c2b = sys.argv[2]
out_dir = sys.argv[3]
# Comma separated string
bin_ignore_list = sys.argv[4].strip('\'').split(',')
header = eval(sys.argv[5])
binner = sys.argv[6]

if binner:
    binner += '_'

assembly_dict = {}

load_fasta(assembly, assembly_dict)

# Open binner output.
with open(binner_c2b, 'r') as b:
    if header:
        next(b)
    # For each line get contig and bin it was put into.
    for line in b:
        line = line.strip('\t\n ')
        line_contig = line.split()[0]
        line_bin = line.split()[1]
        if line_bin not in bin_ignore_list:
            bin_file_name = f'{binner}{line_bin}.fasta'
            out_bin = f'{out_dir}/{bin_file_name}'
            cont_seq = assembly_dict[line_contig]
            with open(out_bin, "a") as out:
                out.write(f'>{line_contig}\n{cont_seq}\n')
