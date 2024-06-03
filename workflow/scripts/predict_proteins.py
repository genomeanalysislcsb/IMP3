import argparse
import gzip
import multiprocessing.pool
import pyrodigal

def run_gene_finder_on_sequence(args):
    finder, sequence_id, sequence = args
    return (sequence_id, finder.find_genes(sequence))

def predict_proteins(input_file, output_file, num_cores, write_gff=False):
    gene_finder = pyrodigal.GeneFinder(meta=True, mask=True, min_mask=0)

    output_file_base = output_file.split("/")[-1]
    genome_id = output_file_base.split("_protein.faa.gz")[0]

    # Read the input genomic FASTA file
    print(f"{genome_id}: Predicting proteins from {input_file} ...")
    with gzip.open(input_file, "rt") as f:
        headers = []
        sequences = []
        current_seq = []
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip('\n >').replace(" ", "_"))
                if current_seq:
                    # Some sequences have LOTS of trailing Ns (???), remove them
                    sequences.append("".join(current_seq).rstrip('N'))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            # Some sequences have LOTS of trailing Ns (???), remove them
            sequences.append("".join(current_seq).rstrip('N'))

    # Run Pyrodigal prediction in parallel
    print(f"{genome_id}: Running Pyrodigal prediction on {len(headers)} sequences with {num_cores} cores ...")
    finder2id2seq =  [(gene_finder, header, sequence) for header, sequence in zip(headers, sequences)]
    with multiprocessing.pool.ThreadPool(num_cores) as pool:
        predictions = pool.map(run_gene_finder_on_sequence, finder2id2seq)

    # Write the predicted proteins to the output file
    print(f"{genome_id}: Writing predicted proteins to {output_file} ...")
    with gzip.open(output_file, "wt") as f:
        for contig_id, genes in predictions:
            genes.write_translations(f, sequence_id=contig_id, include_stop=False)

    if write_gff:
        # Also write gff file
        gff_file = input_file.replace("_genomic.fna.gz", ".gff.gz")
        print(f"{genome_id}: Writing GFF file to {gff_file} ...")
        with gzip.open(gff_file, "wt") as f:
            for contig_id, genes in predictions:
                genes.write_gff(f, sequence_id=contig_id, header=False)

def main():
    parser = argparse.ArgumentParser(description="Predict proteins from genomic sequences using Pyrodigal.")
    parser.add_argument("--input", required=True, help="Input gzipped genomic FASTA file.")
    parser.add_argument("--output", required=True, help="Output gzipped protein FASTA file.")
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use.")
    parser.add_argument("--write_gff", action="store_true", help="Write GFF file.")
    
    args = parser.parse_args()
    
    predict_proteins(args.input, args.output, args.cores, args.write_gff)

if __name__ == "__main__":
    main()
