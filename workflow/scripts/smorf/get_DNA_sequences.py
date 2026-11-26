# smORF module of IMP3. From the smORF prediction pipeline of Luis F. Delgado and Julia Ortis Sunyer
# (gitlab.com/uniluxembourg/lcsb/systems-ecology/smorf_methods_paper; Davin, Ortis Sunyer, Delgado et al.,
# Nat Commun 2026, doi:10.1038/s41467-026-72762-5). GPL-3.0 (LICENSE.GPL-3.0 in this folder).
import os
import argparse
import pandas as pd

def get_proteins_ids(filein):
    # >bhc01::Sbhc01_merged_C0000000022_1 ...
    proteins = dict()
    with open(filein, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                header = line[1:].split()[0]  # bhc01::Sbhc01_merged_C0000000022_1
                if "::" in header:
                    prot = header.split("::")[1]  # Sbhc01_merged_C0000000022_1
                else:
                    prot = header
                proteins[prot] = header
    return proteins  



def get_gff(proteins, filein1, filein2, outfile, extra=()):
    selected_lines = []
    contigs = set("_".join(k.split("_")[:-1]) for k in proteins.keys())

    for filein in [filein1, filein2] + list(extra):
        with open(filein, "r") as fin:
            for line in fin:
                line = line.rstrip()
                if not line.startswith("#"):
                    fields = line.split("\t")
                    if len(fields) != 9:
                        continue  # skip malformed lines

                    if fields[0] in contigs:
                        tipo = fields[0].split("_")[0]
                        attr_parts = fields[8].split(";")

                        for part in attr_parts:
                            if part.startswith("ID="):
                                raw_id = part.split("ID=")[1]
                                new_id = f"{tipo}_{raw_id}"

                                if new_id in proteins:
                                    # Replace only the ID field
                                    updated_parts = [
                                        f"ID={proteins[new_id]}" if p.startswith("ID=") else p
                                        for p in attr_parts
                                    ]
                                    fields[8] = ";".join(updated_parts)
                                    selected_lines.append(fields)
                                    break  # Only need to match once

    # Create DataFrame
    df = pd.DataFrame(selected_lines, columns=["#seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

    # Ensure numeric sorting
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    # After creating df and before sorting:
    df["ORIGINAL"] = df["#seqname"].str.split("_", n=1).str[1]


    # Sort and write output
    df = df.sort_values(by=["ORIGINAL", "start"], ascending=[True, True]).reset_index(drop=True)
    df = df.drop(columns=["ORIGINAL"])

    df.to_csv(outfile, sep="\t", index=False)


def print_DNA_seqs(proteins, filein, fileout):
    # >Sbhc01_merged_C0000000005_1 ...
    copy = False
    header = ""
    seq = ""
    with open(fileout, "w") as fout, open(filein, "r") as fin:  
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"): 
                if copy:
                    print(f"{header}\n{seq}", file=fout)
                    copy = False
                    seq = ""
                    header = ""
                ID = line[1:].split()[0]
                if ID in proteins:
                    header = line.replace(ID, proteins[ID])  
                #    proteins.pop(ID, None) 
                    copy = True
                else:
                    copy = False
            elif copy:
                seq += line
        if copy:
            print(f"{header}\n{seq}", file=fout)
    return proteins        

def print_DNA_seqs2(proteins, filein, fileout):
    # >Sbhc01_merged_C0000000005_1 ...
    copy = False
    header = ""
    seq = ""
    with open(fileout, "a") as fout, open(filein, "r") as fin:  
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"): 
                if copy:
                    print(f"{header}\n{seq}", file=fout)
                    copy = False
                    seq = ""
                    header = ""
                ID = line[1:].split()[0]
                if ID in proteins:
                    header = line.replace(ID, proteins[ID])  
                #    proteins.pop(ID, None) 
                    copy = True
                else:
                    copy = False
            elif copy:
                seq += line
        if copy:
            print(f"{header}\n{seq}", file=fout)

def main():
    parser = argparse.ArgumentParser(
        description='This program extracts the DNA sequences based on the protein ID',
        usage='python get_DNA_sequences.py [options]'
    )

    parser.add_argument('-gd', dest='gd', help='DNA fasta file default (.fna)', required=True)
    parser.add_argument('-gs', dest='gs', help='DNA fasta file small (.fna)', required=True)  
    parser.add_argument('-ad', dest='ad', help='gff file default (.gff)', required=True)
    parser.add_argument('-asm', dest='asm', help='gff file small (.gff)', required=True)  
    parser.add_argument('-p', dest='p', help='Selected proteins .fasta file', required=True)
    parser.add_argument('-o', dest='o', help='Output file for selected genes', required=True)
    parser.add_argument('-a', dest='a', help='Output gff file for selected genes', required=True)
    # IMP3: further gene sets in the same format (e.g. the small CDS of the IMP3 annotation, prefix imp3_)
    parser.add_argument('-gx', dest='gx', nargs='*', default=[], help='further DNA fasta files (.fna)')
    parser.add_argument('-ax', dest='ax', nargs='*', default=[], help='further gff files, same order as -gx')

    args = parser.parse_args()

    # Load potential smORF protein IDs
    proteins_ids = get_proteins_ids(args.p) 
    proteins_ids=print_DNA_seqs(proteins_ids, args.gd, args.o)
    print_DNA_seqs2(proteins_ids, args.gs, args.o)
    for gx in args.gx:
        print_DNA_seqs2(proteins_ids, gx, args.o)
    #print(proteins_ids.keys())
    get_gff(proteins_ids, args.ad, args.asm, args.a, args.ax)


if __name__ == "__main__":
    main()
