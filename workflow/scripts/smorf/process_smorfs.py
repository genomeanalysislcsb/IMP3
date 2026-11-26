# smORF module of IMP3. From the smORF prediction pipeline of Luis F. Delgado and Julia Ortis Sunyer
# (gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper; Davin, Ortis Sunyer, Delgado et al.,
# Nat Commun 2026, doi:10.1038/s41467-026-72762-5). GPL-3.0 (LICENSE.GPL-3.0 in this folder).
import os
import argparse


def process_sequences(file_path, output_file, sample=None, max_len=100):
    sequences = {}
    current_sequence_name = False
    current_sequence = ""

    with open(file_path, "r") as fasta_in, open(output_file, "w") as fasta_out:
        for line in fasta_in:
            line = line.rstrip()
            if line.startswith(">"):
                if current_sequence_name:
                    if len(current_sequence) <= max_len and not "rbs_motif=None" in current_sequence_name and "partial=00" in current_sequence_name:
                        header = f"{current_sequence_name};smORF [len {len(current_sequence)}]"
                        if not ("::") in header:

                            sample = sample or os.path.basename(file_path).split(".")[0]
                            header = f"{sample}::{header}"    
 
                        print(f">{header}\n{current_sequence}", file=fasta_out)

                    
                current_sequence_name = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence_name:
#            if "partial=00" in current_sequence_name and "rbs_motif=None" not in current_sequence_name and len(current_sequence) <= max_len:
            if len(current_sequence) <= max_len and not "rbs_motif=None" in current_sequence_name and "partial=00" in current_sequence_name:
                            header = f"{current_sequence_name};smORF [len {len(current_sequence)}]"
                            if not ("::") in header:
                                sample = sample or os.path.basename(file_path).split(".")[0]
                                header = f"{sample}::{header}" 
                            print(f">{header}\n{current_sequence}", file=fasta_out)


def main():
    parser = argparse.ArgumentParser(description="Select smORFs.")
    parser.add_argument("-i", required=True, help="input fasta file")
    parser.add_argument("-o",required=True, help="file to store the output files.")
    # IMP3: the sample name is given (sample names may contain dots), and the maximum length (aa) is a parameter
    parser.add_argument("-s", default=None, help="sample name (default: input file name up to the first dot)")
    parser.add_argument("-m", type=int, default=100, help="maximum length of a smORF in amino acids (default 100)")

    args = parser.parse_args()

    process_sequences(args.i, args.o, args.s, args.m)

if __name__ == "__main__":
    main()


