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
USER = snakemake.params.user
PASSW = snakemake.params.password
PORT = snakemake.params.port

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


#insert contig dictionaries into MongoDB
MongoClient('mongodb://'+ USER + ':' + PASSW + '@localhost:' + PORT + '/admin')
db = client[LIB]
print("Filling MAG information into the database.")
collection = db.bins
oldCollSize = db.bins.count()
for bin in magdict:
    magdict[bin]['bin'] = bin
    db.bins.insert_one(magdict[bin])
newCollSize = db.bins.count()
print("MAGS inserted into database:", newCollSize - oldCollSize)
print("there are now", newCollSize, "documents in the MAG collection.")
client.close()
