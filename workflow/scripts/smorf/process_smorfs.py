import os
import argparse


def process_sequences(file_path, output_file):
    sequences = {}
    current_sequence_name = False
    current_sequence = ""

    with open(file_path, "r") as fasta_in, open(output_file, "w") as fasta_out:
        for line in fasta_in:
            line = line.rstrip()
            if line.startswith(">"):
                if current_sequence_name:
                    if len(current_sequence) <= 100 and not "rbs_motif=None" in current_sequence_name and "partial=00" in current_sequence_name:
                        header = f"{current_sequence_name};smORF [len {len(current_sequence)}]"
                        if not ("::") in header:

                            sample = os.path.basename(file_path).split(".")[0]
                            header = f"{sample}::{header}"    
 
                        print(f">{header}\n{current_sequence}", file=fasta_out)

                    
                current_sequence_name = line[1:]
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence_name:
#            if "partial=00" in current_sequence_name and "rbs_motif=None" not in current_sequence_name and len(current_sequence) <= 100:
            if len(current_sequence) <= 100 and not "rbs_motif=None" in current_sequence_name and "partial=00" in current_sequence_name:
                            header = f"{current_sequence_name};smORF [len {len(current_sequence)}]"
                            if not ("::") in header:
                                sample = os.path.basename(file_path).split(".")[0]
                                header = f"{sample}::{header}" 
                            print(f">{header}\n{current_sequence}", file=fasta_out)


def main():
    parser = argparse.ArgumentParser(description="Select smORFs.")
    parser.add_argument("-i", required=True, help="input fasta file")
    parser.add_argument("-o",required=True, help="file to store the output files.")

    args = parser.parse_args()

    process_sequences(args.i, args.o)

if __name__ == "__main__":
    main()


