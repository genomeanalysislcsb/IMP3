#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob


def write_bin_table(binner_fasta_dir, binner_output_file):
    with open(binner_output_file, 'w') as of:
        for bin_fasta in glob.glob(binner_fasta_dir + '/*.fasta'):
            bin_name = bin_fasta.replace('.fasta', '').split('/')[-1]
            with open(bin_fasta, 'r') as f:
                for line in f:
                    line = line.strip('\n \t')
                    if line.startswith('>'):
                        contig = line.strip('>')
                        of.write(contig + '\t' + bin_name + '\n')


fasta_dir = sys.argv[1]
output_file = sys.argv[2]

write_bin_table(fasta_dir, output_file)
