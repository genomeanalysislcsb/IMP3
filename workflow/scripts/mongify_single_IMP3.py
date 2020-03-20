#!/usr/bin/python

import os
import sys
import argparse
import re
import numpy
import collections
from pymongo import MongoClient

LIB = snakemake.params.sample
ASS = snakemake.params.ass
BINNERS = snakemake.params.binners.split()
DBS = snakemake.params.dbs.split()
#SDIR = should be current working directory

#function to be used to simplify lists of taxa
def enumerateTaxa(taxList):
    kp = collections.Counter(taxList).most_common()
    init = 1
    taxString = ""
    for tax,occ in kp:
        if init == 1:
            init = 0
            taxString = tax + "(" + str(occ) + ");"
        else:
            taxString += tax + "(" + str(occ) + ");"
    taxString = taxString.rstrip(";")
    return taxString


# read file with selected MAG information
## here we want to have the MAG name, taxonomy, quality and GRiD result + sample of origin and method;
## all MAG will become one collection 
magdict = {}

binFile = "Stats/%s/%s.bins.tsv" % (ASS,ASS)
print("Reading MAG information from ", binFile)
header = 1
bin_file = open(binFile, "r")
while 1:
    line = bin_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0: selected_by_DASTool	1: binner 2: bin 3: uniqueBacSCGs.binner 4: redundantBacSCGs.binner 5: uniqueArcSCGs.binner 6: redundantArcSCGs.binner 7: bacRatio.binner
                               # 8: arcRatio.binner 9: size.binner 10: contigs.binner 11: N50.binner 12: binScore.binner 13: SCG_completeness.binner
                               # 14: SCG_redundancy.binner 15: uniqueBacSCGs.DASTool 16: redundantBacSCGs.DASTool 17: uniqueArcSCGs.DASTool 18: redundantArcSCGs.DASTool 19: bacRatio.DASTool
                               # 20: arcRatio.DASTool 21: size.DASTool 22: contigs.DASTool 23: N50.DASTool 24: binScore.DASTool 25: SCG_completeness.DASTool 26: SCG_redundancy.DASTool
                               # 27: GRiD 28: X95..CI 29: GRiD.unrefined 30: Species.heterogeneity 31: Coverage 32: dnaA.ori.ratio 33: ter.dif.ratio
                               # (GTDBtk) 34: classification 35: fastani_reference 36: fastani_reference_radius 37: fastani_taxonomy 38: fastani_ani 39: fastani_af 40: closest_placement_reference
                               # 41: closest_placement_taxonomy 42: closest_placement_ani 43: closest_placement_af 44: pplacer_taxonomy 45: classification_method 46: note 47: other_related_references.genome_id.species_name.radius.ANI.AF.
                               # 48: aa_percent 49: translation_table 50: red_value 51: warnings
        if tab[0] != "":
            curr_bin = tab[0]
            if len(tab) > 34:
                magdict[curr_bin] = {'sample': LIB, 'MAG': curr_bin, 'binScore': float(tab[24]), 'binningMethod': tab[1],
                                     'GRiD': float(tab[27]) , 'depth': float(tab[31]) , 'dnaAToOri': float(tab[32]) , 'terToDif': float(tab[33]),
                                     'taxstring': tab[34], 'taxMethod': tab[45]}
            else:
                magdict[curr_bin] = {'sample': LIB, 'MAG': curr_bin, 'binScore': float(tab[24]), 'binningMethod': tab[1],
                                     'GRiD': float(tab[27]) , 'depth': float(tab[31]) , 'dnaAToOri': float(tab[32]) , 'terToDif': float(tab[33])}
            
bin_file.close()
print("Read MAG summary information")

# read file of contig information and make dictionary with all contigs
contig_dict = {}
cct = 0
contigFile = "Stats/%s.assembly.length.txt" % ASS
print("Reading contig information from ", contigFile)
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    cct += 1
    line = line.rstrip()
    tab = line.split("\t") # 0:sequenceID, 1:length
    contig_dict[tab[0]] = {'sample': LIB,'length' : int(tab[1]), 'GCperc' : 0, 'aveCov' : 0 }
contig_file.close()
print("Read contig length information on ", cct, " contigs.")

# read file of contig GC content and update dictionary 
header = 1
cct = 0
contigFile = "Stats/%s.assembly.gc_content.txt" % ASS
print("Reading contig information from ", contigFile)
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        cct += 1
        line = line.rstrip()
        tab = line.split("\t") # 0:sequenceID, 1:GC
        contig_dict[tab[0]]['GCperc'] = float(tab[1])
contig_file.close()
print("Read contig GC information on ", cct, " contigs.")


# read file of DNA coverage of the same sample for contigs and update contig dictionary
covDNAcontigFile = "Stats/%s/mg.assembly.contig_depth.txt" % ASS
if os.path.isdir(covDNAcontigFile):
    cct = 0
    print("Reading coverage data for DNA on contigs from ", covDNAcontigFile)
    covDNAcontig_file = open(covDNAcontigFile, "r")
    while 1:
        line = covDNAcontig_file.readline()
        if line == "":
            break
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Average.coverage 
        contig_dict[tab[0]]['aveCov'] = float(tab[1])
    covDNAcontig_file.close()
    print("Read contig MG coverage depth information on ", cct, " contigs.")
    
print("gathered basic information on ", len(contig_dict), "contigs.")

# read coordinates
coordFile = "Binning/%s.vizbin.with-contig-names.points" % ASS
if os.path.isdir(coordFile):
    vizCnt = 0
    print("Reading contig BHSNE coordinates from ", coordFile)
    coord_file = open(coordFile, "r")
    while 1:
        line = coord_file.readline()
        if line == "":
            break
        else:
            vizCnt += 1
            line = line.rstrip()
            tab = line.split("\t") # 0:contig 1:x, 2:y
            contig_dict[tab[0]]['coords'] = [float(tab[1]),float(tab[2])] 
    coord_file.close()
    print("gathered vizbin coordinates of ", vizCnt, "contigs")

#read single binner memberships
for binner in BINNERS:
    clusterFile = "/Binning/%s/scaffold2bin.tsv" % binner
    if os.path.isdir(clusterFile):
        contigCnt = 0
        clusterName = binner + "_bin"
        print("Reading contig cluster membership from ", clusterFile)
        cluster_file = open(clusterFile, "r")
        while 1:
            line = cluster_file.readline()
            if line == "":
                break
            contigCnt +=  1
            line = line.rstrip()
            tab = line.split("\t") # 0:contig, 1:class (contig membership)
            contig_dict[tab[0]][clusterName] = tab[1]
        cluster_file.close()
        print("gathered ", binner, " information on ", contigCnt, "contigs")


#read DAStool
dastoolFile = "Binning/selected_DASTool_scaffolds2bin.txt"
if os.path.isdir(dastoolFile):
    dastoolCnt = 0
    print("Reading contig bin membership from ", dastoolFile)
    dastool_file = open(dastoolFile, "r")
    while 1:
        line = dastool_file.readline()
        if line == "":
            break
        dastoolCnt +=  1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['DASTool_MAG'] = tab[1]
    dastool_file.close()
    print("read DAStool information on ", dastoolCnt, "contigs")


# read Kraken data and update contig dictionary
krakenFile = "Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" % ASS
if os.path.isdir(krakenFile):
    header = 1
    krakCnt = 0
    print("Reading contig Kraken information from ", krakenFile)
    kraken_file = open(krakenFile, "r")
    while 1:
        line = kraken_file.readline()
        if line == "":
            break
        line = line.rstrip()
        tab = line.split("\t") # 0:contig 1: length 2: entropy 3: annotationLevel 4+: rest is free
        if header == 1:
            header = 0
            levels = ["kraken_" + s for s in tab]
        else:
            krakCnt += 1
            contig_dict[tab[0]]['krakenAnnotationLevel'] = tab[3]
            contig_dict[tab[0]].update(dict(filter(lambda x:x[1]!="NA", zip(levels,tab))))
    kraken_file.close()
    print("read kraken information on ", krakCnt, "contigs")

protein_dict = {}
rrna_dict = {}
trna_dict = {}
# read gff file
protFile = "Analysis/annotation/annotation_CDS_RNA_hmms.gff" if os.path.isdir("Analysis/annotation/annotation_CDS_RNA_hmms.gff") else "Analysis/annotation/annotation.filt.gff"
if os.path.isdir(protFile):
    pct = 0
    rct = 0
    tct = 0
    lct = 0
    print("Reading protein information from ", protFile)
    prot_file = open(protFile, "r")
    while 1:
        line = prot_file.readline()
        if line == "":
            break
        else:
            lct += 1
            line = line.rstrip()
            tab = line.split("\t") # 0:contig, 1:source, 2:kind, 3:start, 4:end, 5:.n, 6:sense, 7:0, 8:attributes
            if len(tab) >= 9:
                atts = tab[8].split(";")
                for att in atts:
                    if att.startswith("ID="):
                        gene = att.replace("ID=","")
                    if tab[2] == "CDS":
                        if att.startswith("partial="):
                            part = att.replace("partial=","")
                            completeness = "complete" if part == "00" else "incomplete"
                    elif tab[2] == "rRNA":
                        if att.startswith("product="):
                            completeness = "incomplete" if att.endswith(" (partial)") else "complete"
                            prod_t = att.replace("product=","")
                            rRNA_type = prod_t.replace(" ribosomal RNA (partial)","") if att.endswith(" (partial)") else att.replace(" ribosomal RNA","")
                    elif tab[2] == "tRNA" or tab[2] == "tmRNA":
                        if att.startswith("product="):
                            tRNA_type = att.replace("product=","")
                if tab[3] != "" and tab[4]!= "":
                    if tab[2] == "CDS":
                        pct += 1
                        protein_dict[gene] = {'contig' : tab[0], 'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]),
                                              'completeness': completeness,'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0, 'readsDNA' : 0, 'readsRNA': 0}
                        protein_dict[gene].update(dict([item for sublist in [[(db[:-1],att.replace(db,"")) for db in [s + "=" for s in dbs] if att.startswith(db)] for att in atts] for item in sublist]))
                    elif tab[2] == "rRNA":
                        rct += 1
                        rrna_dict[gene] = {'contig' : tab[0], 'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]),
                                           'kind' : tab[2], 'rRNA': rRNA_type,'aveCovDNA' : 0, 'aveCovRNA' : 0}
                    elif tab[2] == "tRNA" or "tmRNA":    
                        tct += 1
                        trna_dict[gene] = {'contig' : tab[0], 'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]),
                                           'kind' : tab[2], 'tRNA': tRNA_type, 'aveCovDNA' : 0, 'aveCovRNA' : 0}
                else:
                    print("Missing genomic coordinates for ", gene)
            else:
                print("missing attributes for the gene in line ", lct, "of ",  protFile)
    prot_file.close()
    print("read information on ", lct, " regions: ", pct, " proteins, ", rct, " rRNAs, ", tct, " tRNAs/tmRNAs.")

# read file of DNA coverage for genes and update gene dictionary
covDNAgeneFile = "Stats/mg/annotation/mg.gene_depth.avg"
if os.path.isdir(covDNAgeneFile):
    gct = 0
    print("Reading coverage data for DNA on genes from ", covDNAgeneFile)
    covDNAgene_file = open(covDNAgeneFile, "r")
    while 1:
        line = covDNAgene_file.readline()
        if line == "":
            break
        gct += 1
        line = line.rstrip()
        tab = line.split("\t") # 0:geneID,	1:average depth
        if protein_dict.has_key(tab[0]):
            protein_dict[tab[0]]['aveCovDNA'] = float(tab[1])
        elif rrna_dict.has_key(tab[0]):
            rrna_dict[tab[0]]['aveCovDNA'] = float(tab[1])
        elif trna_dict.has_key(tab[0]):
            trna_dict[tab[0]]['aveCovDNA'] = float(tab[1])
    covDNAgene_file.close()
    print("read DNA coverage information for ", gct, " genes.")

# read file of RNA coverage for genes and update gene dictionary
covRNAgeneFile = "Stats/mt/annotation/mt.gene_depth.avg"
if os.path.isdir(covRNAgeneFile):
    gct = 0
    print("Reading forward coverage data for RNA on genes from ", covRNAgeneFile)
    covRNAgene_file = open(covRNAgeneFile, "r")
    while 1:
        line = covRNAgene_file.readline()
        if line == "":
            break
        gct += 1
        line = line.rstrip()
        tab = line.split("\t") # 0:geneID,	1:average depth
        if protein_dict.has_key(tab[0]):
            protein_dict[tab[0]]['aveCovRNA'] = float(tab[1])
        elif rrna_dict.has_key(tab[0]):
            rrna_dict[tab[0]]['aveCovRNA'] = float(tab[1])
        elif trna_dict.has_key(tab[0]):
            trna_dict[tab[0]]['aveCovRNA'] = float(tab[1])
    covRNAgene_file.close()
    print("read RNA coverage information for ", gct, " genes.")

# read files of RNA reads per gene and update protein dictionary
readsDNAgeneFile = "Analysis/annotation/mg.CDS_counts.tsv"
if os.path.isdir(readsDNAgeneFile):
    gct = 0
    print("Reading coverage data for DNA on genes from ", covDNAgeneFile)
    header = 2
    readsDNAgene_file = open(readsDNAgeneFile, "r")
    while 1:
        line = readsDNAgene_file.readline()
        if line == "":
            break
        if header > 0:
            header -= 1
        else:
            gct += 1
            line = line.rstrip()
            tab = line.split("\t") # 0:Geneid  1:Chr 2:Start  3:End  4:Strand  5:Length 6:Assembly/mg.reads.sorted.ba
            if protein_dict.has_key(tab[0]):
                protein_dict[tab[0]]['readsDNA'] = float(tab[6])
    readsDNAgene_file.close()
    print("read DNA coverage information for ", gct, " genes.")

# read files of RNA reads per gene and update protein dictionary
readsRNAgeneFile = SDIR + "/Analysis/annotation/mt.CDS_counts.tsv"
if os.path.isdir(readsRNAgeneFile):
    print("Reading forward coverage data for RNA on genes from ", covRNAgeneFile)
    header = 2
    readsRNAgene_file = open(readsRNAgeneFile, "r")
    while 1:
        line = readsRNAgene_file.readline()
        if line == "":
            break
        if header > 0:
            header -= 1
        else:
            gct += 1
            line = line.rstrip()
            tab = line.split("\t") # 0:Geneid  1:Chr 2:Start  3:End  4:Strand  5:Length 6:Assembly/mt.reads.sorted.bam
            if protein_dict.has_key(tab[0]):
                protein_dict[tab[0]]['readsRNA'] = float(tab[6])
    readsRNAgene_file.close()
    print("read RNA coverage information for ", gct, " genes.")


#insert protein information into the contig dictionary
contcount=0
if protein_dict:
    print("inserting information on ", len(protein_dict), "genes into contig data")
    for prot in protein_dict:
        protein_dict[prot]['gene'] = prot
        contig = protein_dict[prot]['contig']
        if contig not in contig_dict:
            print("no contig information available for protein", prot, " on ", contig)
        else:
            if "genes" not in contig_dict[contig]:
                contig_dict[contig]['genes'] = [protein_dict[prot]]
                contcount += 1
            else:
                contig_dict[contig]['genes'].append(protein_dict[prot])
    protein_dict = {}

if rrna_dict:
    print("inserting information on ", len(rrna_dict), "rRNA genes into contig data")
    for rna in rrna_dict:
        rrna_dict[rna]['gene'] = rna
        contig = rrna_dict[rna]['contig']
        if contig not in contig_dict:
            print("no contig information available for rRNA", rna, " on ", contig)
        else:
            if "genes" not in contig_dict[contig]:
                contig_dict[contig]['genes'] = [rrna_dict[rna]]
                contcount += 1
            else:
                contig_dict[contig]['genes'].append(rrna_dict[rna])
    rrna_dict = {}

if trna_dict:
    print("inserting information on ", len(trna_dict), "rRNA genes into contig data")
    for rna in trna_dict:
        trna_dict[rna]['gene'] = rna
        contig = trna_dict[rna]['contig']
        if contig not in contig_dict:
            print "no contig information available for tRNA", rna, " on ", contig
        else:
            if "genes" not in contig_dict[contig]:
                contig_dict[contig]['genes'] = [trna_dict[rna]]
                contcount += 1
            else:
                contig_dict[contig]['genes'].append(trna_dict[rna])
    trna_dict = {}

if contcount > 0:
    print("gene information inserted for ", contcount, " contigs")
else:
    print("No gene information.")


#insert contig dictionaries into MongoDB
client = MongoClient()
db = client[LIB]
print("Filling contig information into the database.")
oldCollSize = db.contigs.count()
for cont in contig_dict:
    contig_dict[cont]['contig'] = cont
    db.contigs.insert_one(contig_dict[cont])
newCollSize = db.contigs.count()
print("contigs inserted into database:", newCollSize - oldCollSize)
print("there are now", newCollSize, "documents in the contig collection.")
print("Filling MAG information into the database.")
oldCollSize = db.bins.count()
for bin in magdict:
    magdict[bin]['bin'] = bin
    db.bins.insert_one(magdict[bin])
newCollSize = db.bins.count()
print("MAGS inserted into database:", newCollSize - oldCollSize)
print("there are now", newCollSize, "documents in the MAG collection.")
