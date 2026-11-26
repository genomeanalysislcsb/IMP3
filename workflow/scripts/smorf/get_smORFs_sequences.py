import os
import argparse
from collections import defaultdict
from itertools import combinations


def compare_sources(sources: dict):
    """
    Compare all pairs of sets in `sources` dictionary and print shared and unique item counts.
    """
    for (k1, s1), (k2, s2) in combinations(sources.items(), 2):
        shared = s1 & s2
        only_in_s1 = s1 - s2
        only_in_s2 = s2 - s1
        print(f"{k1} ∩ {k2}: {len(shared)} shared items")
        print(f"{k1} only: {len(only_in_s1)} unique items")
        print(f"{k2} only: {len(only_in_s2)} unique items")
        print("-" * 50)


def get_aragorn(fasta_path, data):
    valid_ids = set()

    with open(fasta_path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Look for header lines starting with ">" and containing "::"
        if line.startswith(">") and "::" in line:
            full_header = line
            i += 1
            # Check if the next line is NOT "0 genes found"
            if i < len(lines) and "0 genes found" not in lines[i]:
                # Extract only the part after ">" up to the first space
                seq_id = full_header[1:].split(" ")[0]
                valid_ids.add(seq_id)
        else:
            i += 1


    print(f"INFO: Total spurius smORF - {fasta_path} : {len(valid_ids)}") 
    data.update(valid_ids)

    return data, valid_ids


def get_tRNAscan(file_path, data):
    counter = set()

    with open(file_path, 'r') as f:
        lines = f.read().splitlines()[3:]  # Skip headers and dashes

    for line in lines:
        parts = line.strip().split()
        if parts:
            counter.add(parts[0])

    print(f"INFO: Total spurious smORF - {file_path} : {len(counter)}") 
    
    data.update(counter)
    return data, counter


def get_barrnap(file_path, data):
    counter = set()
    """
    ##gff-version 3
    bhc01::Sbhc01_merged_C0000045036_1	barrnap:0.9	rRNA	10	105	5e-08	+	.	Name=5S_rRNA;product=5S ribosomal RNA
    """
    with open(file_path, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if not line.startswith("#"):
                id=line.split("\t")[0]
                counter.add(id)

    print(f"INFO: Total spurius smORF - {file_path} : {len(counter)}") 
    
    data.update(counter)
    return data, counter

def get_RFAM_hits(file_path, data):
    counter = set()
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) > 3:
                counter.add(parts[3]) 
    print(f"INFO: Total spurius smORF - {file_path} : {len(counter)}") 
    
    data.update(counter)           
    return data, counter

def get_antiFAM_hits(file_path, data):
    counter=set()
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) > 3:
                counter.add(parts[0]) 

    print(f"INFO: Total spurius smORF - {file_path} : {len(counter)}") 
    
    data.update(counter)            

    return data, counter


def print_gff(toremove, filein, fileout):
    with open(filein, "r") as fin, open(fileout, "w") as fout:
        for line in fin:
            line = line.rstrip()
            if line.startswith("#"):
                print(line, file = fout)
            else:
                fields = line.split("\t")
                id=fields[8].split(";")[0]    
                id=id.replace("ID=","")
                if not id in toremove:
                    print(line, file = fout)



def print_seqs(toremove, filein, fileout, log):
    counterI=0
    counterF=0
    copy = False
    header = ""
    seq = ""
    with open(fileout, "w") as fout, open(filein, "r") as fin:  
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"): 
                counterI += 1
                if copy:
                    counterF += 1
                    print(f"{header}\n{seq.replace('*','')}", file=fout)
                ID = line[1:].split()[0]
                if ID not in toremove:
                    header = line
                    copy = True
                    seq = ""
                else:
                    copy = False
            elif copy:
                seq += line
        if copy:
            counterF += 1
            print(f"{header}\n{seq.replace('*','')}", file=fout)
    if log:
        print(f"INFO: Total number of initial sequences: {counterI}")
        print(f"INFO: Total number of final sequences: {counterF}")
        print(f"INFO: Percentage of spurius sequences: {round((counterI-counterF)*100/counterI,2)}%")

def main():
    parser = argparse.ArgumentParser(
        description='This program extracts the DNA sequences based on the protein ID',
        usage='python get_DNA_sequences.py [options]'
    )

    parser.add_argument('-g', help='DNA fasta file (.fna)', required=True)  
    parser.add_argument('-p', help='Selected proteins .fasta file', required=True)
    parser.add_argument('-a', help='Selected proteins .gff file', required=True)
    parser.add_argument('-rfam', help='RFAM result .tbl file', default=None)
    parser.add_argument('-antifam', help='AntiFAM result .tbl file', default=None)
    parser.add_argument('-arbac', help='Aragorn bacteria result file', default=None)
    parser.add_argument('-arstd', help='Aragon std result file', default=None)
    parser.add_argument('-barrnap_bac', help='barrnap_bac result .tbl file', default=None)
    parser.add_argument('-barrnap_arch', help='barrnap_arch result .tbl file', default=None)
    parser.add_argument('-barrnap_mito', help='barrnap_mito result .tbl file', default=None) 
    parser.add_argument('-tRNAscan_arch',  help='tRNAscan_arch result .tbl file', default=None)  
    parser.add_argument('-tRNAscan_bact', help='tRNAscan_bact result .tbl file', default=None)   
    parser.add_argument('-gene', help='Output file for selected genes', required=True)
    parser.add_argument('-prot', help='Output file for selected proteins', required=True)
    parser.add_argument('-gff', help='Output gff file for selected proteins', required=True)
    parser.add_argument('-report', help='Report comparing RFAM vs other sources', default=None)


    args = parser.parse_args()

    # Initialize all sets to None
    rfams = antifams = arbacs = arstds = barrnap_bacs = barrnap_archs = barrnap_mitos = tRNAscan_arch = tRNAscan_bact = None

    to_remove = set()

    # Populate sets from args
    if args.rfam:
        to_remove, rfams = get_RFAM_hits(args.rfam, to_remove)
    if args.antifam:
        to_remove, antifams = get_antiFAM_hits(args.antifam, to_remove)
    if args.arbac:
        to_remove, arbacs = get_aragorn(args.arbac, to_remove)
    if args.arstd:
        to_remove, arstds = get_aragorn(args.arstd, to_remove)
    if args.barrnap_bac:
        to_remove, barrnap_bacs = get_barrnap(args.barrnap_bac, to_remove)
    if args.barrnap_arch:
        to_remove, barrnap_archs = get_barrnap(args.barrnap_arch, to_remove)
    if args.barrnap_mito:
        to_remove, barrnap_mitos = get_barrnap(args.barrnap_mito, to_remove)
    if args.tRNAscan_arch:
        to_remove, tRNAscan_arch = get_tRNAscan(args.tRNAscan_arch, to_remove)
    if args.tRNAscan_bact:
        to_remove, tRNAscan_bact = get_tRNAscan(args.tRNAscan_bact, to_remove)

    print(f"INFO: Total spurious smORF: {len(to_remove)}")


    sources = defaultdict(set)
    if rfams: sources["rfams"] = rfams
    if antifams: sources["antifams"] = antifams

    # Add combined groups
#    sources["Aragorns"] = (arbacs or set()) | (arstds or set())
#    sources["barrnaps"] = (barrnap_bacs or set()) | (barrnap_archs or set()) | (barrnap_mitos or set())
#    sources["tRNAscan"] = (tRNAscan_arch or set()) | (tRNAscan_bact or set())
    sources["no_rfam"] = (barrnap_bacs or set()) | (barrnap_archs or set()) | (barrnap_mitos or set()) | (arbacs or set()) | (arstds or set()) | (antifams or set()) | (tRNAscan_arch or set()) | (tRNAscan_bact or set())

    compare_sources(sources)

    print_seqs(to_remove, args.g, args.gene, log=False)
    print_seqs(to_remove, args.p, args.prot, log=True)
    print_gff(to_remove, args.a, args.gff)

if __name__ == "__main__":
    main()
