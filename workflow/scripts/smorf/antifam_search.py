import argparse
import pyhmmer
import os
from pyhmmer.easel import SequenceFile

#def search_hmm_with_sequence(hmm_file_path, sequence_file_path, user_defined_threshold, num_cpus,output_file):
def search_hmm_with_sequence(hmm_file_path, sequence_file_path, num_cpus,output_file):
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm, SequenceFile(sequence_file_path, digital=True) as sequences, open(output_file, 'wb') as f : 

        targets = sequences.read_block()
    
        results = pyhmmer.hmmsearch(hmm, targets, cpus=num_cpus, bit_cutoffs="gathering") #T=user_defined_threshold
        
        first=True  
        for hits in results:
            if hits.reported:
                if first:
                    hits.write(f, format="targets", header=True)
                    first = False    
                else:    
                    hits.write(f, format="targets", header=False)    

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Parallel HMMER search using pyhmmer")
    
    # Add arguments for HMM files, sequence files, number of CPUs, bit score threshold, and output file
    parser.add_argument("-hmm", "--hmm_file", required=True,
                        help="Paths to the HMM file")
    parser.add_argument("-seqs", "--sequence_file", required=True,
                        help="Paths to the sequence files")
    parser.add_argument("-n", "--num_cpus", type=int, default=os.cpu_count(),
                        help="Number of CPUs to use (default: all available CPUs)")
    #parser.add_argument("-t", "--bit_score_threshold", type=float, required=True,
    #                    help="Bit score threshold for filtering hits")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Output file to write results to")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    #search_hmm_with_sequence(args.hmm_file, args.sequence_file, args.bit_score_threshold, args.num_cpus, args.output_file)
    search_hmm_with_sequence(args.hmm_file, args.sequence_file, args.num_cpus, args.output_file)

if __name__ == "__main__":
    main()
