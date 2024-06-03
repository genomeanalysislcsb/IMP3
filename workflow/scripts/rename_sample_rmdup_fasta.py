#!/bin/python3

#################
### Libraries ###
#################

import argparse

#################
### Functions ###
#################

def read_fasta(fasta_file):
    '''
    Read fasta file and return a dictionary
    Dict is sorted by the length of the sequences
    '''
    fasta = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:]
                fasta[header] = ''
            else:
                fasta[header] += line.strip()

    # Sort the dictionary by the length of its values
    fasta = {k: v for k, v in sorted(fasta.items(), key=lambda item: len(item[1]), reverse=True)}
    return fasta


def read_seqkit_duplicate_file(dup_file):
    '''
    Read seqkit duplicate file and return a dictionary
    The file is tab-separated with the following format:
    'number of duplicates'    'duplicate name1, duplicate name2, ...'
    '''
    dup = {}
    with open(dup_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            representative = line[1].split(', ')[0]
            duplicates = line[1].split(', ')
            dup[representative] = duplicates
    return dup


def rename_fasta(fasta, dup, new_fasta_file, new_dup_file, sample):
    counter = 0
    with open(new_fasta_file, 'w') as ff:
        with open(new_dup_file, 'w') as fd:
            dup_header = 'id\tn_duplicates\tprotein_ids\n'
            fd.write(dup_header)
            for id, seq in fasta.items():
                counter += 1

                new_id = f'{sample}_{counter:09d}'

                # Write to duplicates file
                if id in dup:
                    n_duplicates = len(dup[id])
                    duplicates = ';'.join(dup[id])
                    fd.write(f'{new_id}\t{n_duplicates}\t{duplicates}\n')

                    new_desc = f'n_identical_proteins={n_duplicates}'
                else:
                    fd.write(f'{new_id}\t0\t{id}\n')

                    new_desc = id.split('___')[-1]
                
                # Split sequences longer than 60 characters
                if len(seq) > 60:
                    n = len(seq) // 60
                    seq = '\n'.join([seq[i:i+60] for i in range(0, len(seq), 60)])

                # Write to fasta file
                ff.write(f'>generic|{new_id}|{new_desc}\n{seq}\n')

#################
### Execution ###
#################
                    
# Example: rename_sample_rmdup_fasta.py data/sample.fasta data/sample_duplicates.txt data/sample_rmdup.fasta data/sample_rmdup_duplicates.txt
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename fasta and duplicate files')
    parser.add_argument('sample', help='Sample name')
    parser.add_argument('fasta', help='Input fasta file')
    parser.add_argument('duplicates', help='Input duplicates file')
    parser.add_argument('new_fasta', help='Output fasta file')
    parser.add_argument('new_duplicates', help='Output duplicates file')
    args = parser.parse_args()

    fasta = read_fasta(args.fasta)
    dup = read_seqkit_duplicate_file(args.duplicates)
    rename_fasta(fasta, dup, args.new_fasta, args.new_duplicates, args.sample)


