# smORF module of IMP3. From the smORF prediction pipeline of Luis F. Delgado and Julia Ortis Sunyer
# (gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper; Davin, Ortis Sunyer, Delgado et al.,
# Nat Commun 2026, doi:10.1038/s41467-026-72762-5). GPL-3.0 (LICENSE.GPL-3.0 in this folder).
import os
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import gzip

def identify_compression(file_path):
    magic_dict = {
        b"\x1f\x8b": "gzip",
        b"\x42\x5a": "bz2",
        b"\x50\x4b": "zip",  # PK: ZIP file
        b"\x1f\x9d": "compress (.Z)",
    }

    with open(file_path, "rb") as f:
        file_start = f.read(2)  # Read the first two bytes

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return True

    return False

def process_sample(sample, ids_in_samples, directory_genes, ext):
#def process_sample(sample, ids_in_samples, directory_genes, ext, typein):
    """Process a single sample and return the output."""
    s = sample.split(ext, 1)[0]

    selected_ids = ids_in_samples[s]

    filein = os.path.join(directory_genes, sample)
    result = []

    if identify_compression(filein):
#    if typein == "contig": #and filein.endswith(".gz"):
        fin = gzip.open(filein, "rt")  # Open gzipped file in text mode
    else:
        fin = open(filein, 'r')  # Open regular text file

    with fin:
            write_flag = False
            seq_lines = []
            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    #>bhc01::Sbhc01_merged_C0000000022_1 # 160 # 360 # -1 # ID=21_1;partial=00;start_type=ATG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;gc_cont=0.667
                    if write_flag:
                        result.append(f'{header}\n{"".join(seq_lines)}\n')
                        seq_lines.clear()
                        write_flag = False

                    protid = line.split(None, 1)[0][1:]
                                            
                    if protid in selected_ids:
                        header = line
                        write_flag = True
                elif write_flag:
                    seq_lines.append(line)

            if write_flag:
                result.append(f'{header}\n{"".join(seq_lines)}\n')

    return result

def print_out_dna_sequences(selected_samples, ids_in_samples, directory_genes, output_fasta, num_cpus, extention):
#def print_out_dna_sequences(selected_samples, ids_in_samples, directory_genes, output_fasta, num_cpus, extention, typein):
    """Write DNA sequences to the output file using multiple CPUs."""
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        # Submit tasks for parallel processing
 #           executor.submit(process_sample, sample, ids_in_samples, directory_genes, extention, typein)

        tasks = [
            executor.submit(process_sample, sample, ids_in_samples, directory_genes, extention)
            for sample in selected_samples
        ]

        # Collect and write results
        with open(output_fasta, 'w') as outfile:
            for task in tasks:
                result = task.result()
                outfile.writelines(result)

def get_samples_ids(filein):
    samples = defaultdict(set)
    with open(filein, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                #>bhc01::Sbhc01_merged_C0000312482_1 # 440 # 643 # -1 # ID=310073_1;partial=00;start_type=GTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.520;smORF [len 68] 
                id = line[1:].split(None, 1)[0]
                s = id.split("::", 1)[0]
                samples[s].add(id)
    return samples

def main():
    parser = argparse.ArgumentParser(
        description='This program extracts the DNA sequences based on the protein ID',
        usage='python fromAAtoDNAseqs.py [options]'
    )

    parser.add_argument('-g', dest='g', help='directory with fna fasta files', required=True)
    parser.add_argument('-e', dest='e', default=".fna", help="input fna file extention (defaul: .fna)")
    parser.add_argument('-p', dest='p', help='selected proteins .fasta', required=True)
    parser.add_argument('-o', dest='o', help='output file, selected genes', required=True)
    parser.add_argument('-n', dest='n', type=int, default=1, help='number of CPUs to use (default: 1)')

    args = parser.parse_args()

    # Load selected sample IDs
    selected_samples_ids = get_samples_ids(args.p)
    # Filter relevant samples
    selected_sample_keys = set(selected_samples_ids.keys())


    selected_samples = [
        s for s in os.listdir(args.g)
        if s.endswith(args.e) and s.split(args.e, 1)[0] in selected_sample_keys
    ]

    # Extract sequences using multiple CPUs

    directory = os.path.dirname(args.o)
    
    # Check if the directory exists
    if not os.path.exists(directory):
        # Create the directory if it does not exist
        os.makedirs(directory)            
#    print_out_dna_sequences(selected_samples, selected_samples_ids, args.g, args.o, args.n, args.e, args.t)
    print_out_dna_sequences(selected_samples, selected_samples_ids, args.g, args.o, args.n, args.e)


if __name__ == "__main__":
    main()
