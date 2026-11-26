import re
import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_info(line):
    ID = get_expobiome_id(line)
    s = get_sample(line)
    c = get_contig(line)
    smorf = get_smORF(line)
    source = get_source(line)
    gene = get_gene_id(line)
    return ID, s, c, smorf, source, gene

def get_source(line):
    return line.split("|")[0]

def get_gene_id(line):
    match = re.search(r'gene=([^|]+)', line)
    return match.group(1) if match else "NA"

def get_expobiome_id(line):
    match = re.search(r'ID=([^|]+)', line)
    return match.group(1) if match else "NA"

def get_sample(line):
    match = re.search(r'sample=([^|]+)', line)
    return match.group(1) if match else "NA"

def get_contig(line):
    match = re.search(r'contig=([^|]+)', line)
    return match.group(1) if match else "NA"

def get_smORF(line):
    match = re.search(r'smORF=([^|]+)', line)
    value = match.group(1) if match else "NA"
    return value.replace("__", "")

def read_cluster(filein, fileout):
    from sys import intern  # Save memory by reusing identical strings

    clusters = {}  # {ExpoBiomeID: [cluster_number, representative_ID]}
    current_id = ""
    cluster_number = 0

    header = "Cluster_number\tcluster_membership\tExpoBiomeID\tSource\tSample\tContig\tGene\tIs_smORF"

    with open(filein, "r") as fin, open(fileout, "w") as fout:
        fout.write(header + "\n")
        for line in fin:
            rep, mem = line.rstrip().split("\t")

            # Avoid storing rep metadata unless needed
            ID_rep = get_expobiome_id(rep)
            ID_mem = get_expobiome_id(mem)

            if ID_rep != current_id:
                current_id = ID_rep
                cluster_number += 1
                cluster_label = f"Cluster_{cluster_number}"

            membership = "Representative" if rep == mem else "member"

            # Extract only mem info for output and cluster map
            source = intern(get_source(mem))
            sample = intern(get_sample(mem))
            contig = intern(get_contig(mem))
            gene = intern(get_gene_id(mem))
            smorf = intern(get_smORF(mem))

            clusters[ID_mem] = (cluster_number, ID_rep)

            fout.write(
                f"{cluster_label}\t{membership}\t{ID_mem}\t{source}\t{sample}\t{contig}\t{gene}\t{smorf}\n"
            )

    return clusters

def which_item(line):
    parts = line[1:].split("|")
    for p in range(len(parts)):
        if parts[p].startswith("ID="):
            return p
    return None

def modify_fasta_headers(clusters, filein, dirout):
    fileout = os.path.join(dirout, os.path.basename(filein))
    copy = False
    with open(filein, "r") as fin, open(fileout, "w") as fout:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                try:
                    ID = get_expobiome_id(line[1:])
                except AttributeError:
                    print(f"Could not find ID in line: {line}")
                    continue
                if ID in clusters:
                    i = which_item(line)
                    parts = line[1:].split("|")
                    if i is not None:
                        newID = f"Cluster={clusters[ID][0]}|ID={clusters[ID][1]}"
                        parts[i] = newID
                        new_header = "|".join(parts)
                        new_header = new_header.replace("__","")
                        print(f">{new_header}", file=fout)
                        copy = True
                    else:
                        copy = False
                else:
                    copy = False
            elif copy:
                print(line, file=fout)

def modify_fasta_headers_wrapper(args):
    clusters, filein, dirout = args
    modify_fasta_headers(clusters, filein, dirout)

def main():
    parser = argparse.ArgumentParser(description="Assign identical ExpoBiomeIDs to clustered sequences across multiple samples.")
    parser.add_argument("-tsv", required=True, help="MMseq2 cluster file (.tsv)")
    parser.add_argument("-i", required=True, help="Directory with input FASTA files")
    parser.add_argument("-e", required=True, help="FASTA file extension (e.g., .faa, .fasta)")
    parser.add_argument("-o", required=True, help="Output TSV file with cluster info")
    parser.add_argument("-d", required=True, help="Directory to write modified FASTA files")
    parser.add_argument("-t", type=int, default=4, help="Number of CPUs to use (default: 4)")
    args = parser.parse_args()

    # Create output directories if needed
    if os.path.dirname(args.o):
        os.makedirs(os.path.dirname(args.o), exist_ok=True)
    os.makedirs(args.d, exist_ok=True)

    # Parse clusters
    print(f"INFO: reading {args.tsv}")
    clusters = read_cluster(args.tsv, args.o)

    # Collect FASTA files
    fasta_files = [os.path.join(args.i, f) for f in os.listdir(args.i) if f.endswith(args.e)]
    task_args = [(clusters, f, args.d) for f in fasta_files]

    print(f"INFO: Total number of fasta files to be modified {len(fasta_files)}")
    
    # Run in parallel
    with ProcessPoolExecutor(max_workers=args.t) as executor:
        futures = [executor.submit(modify_fasta_headers_wrapper, a) for a in task_args]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f"Error processing file: {exc}")

if __name__ == "__main__":
    main()
