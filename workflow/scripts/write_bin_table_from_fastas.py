#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import os

def write_bin_table(binner_fasta_dir, binner_output_file):
    with open(binner_output_file, 'w') as of:
        # Check for bin fastas (files with extensions: .fasta, .fa, .fna)
        fasta_extensions = ['*.fasta', '*.fa', '*.fna']
        fasta_files = []
        for ext in fasta_extensions:
            fasta_files += glob.glob(os.path.join(binner_fasta_dir, ext))

        # If no bin fastas are found, exit function since the file is already created empty
        if len(fasta_files) == 0:
            return
        
        # Write the bin table
        for bin_fasta in fasta_files:
            bin_name = os.path.basename(os.path.splitext(bin_fasta)[0])
            try:
                with open(bin_fasta, 'r') as f:
                    for line in f:
                        line = line.strip('\n \t')
                        if line.startswith('>'):
                            contig = line.strip('>')
                            of.write(f"{contig}\t{bin_name}\n")
            except IOError as e:
                print(f"Error reading file {bin_fasta}: {e}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: write_bin_table_from_fastas.py <fasta_dir> <output_file>')
        sys.exit(1)

    fasta_dir = sys.argv[1]
    output_file = sys.argv[2]

    write_bin_table(fasta_dir, output_file)
