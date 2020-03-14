#!/usr/bin/env python3

import os
import sys
import argparse
import operator
import math
import numpy as np

parser = argparse.ArgumentParser(description='Process output of Kraken run on contigs.')
parser.add_argument('inputFile', help='output from Kraken')
parser.add_argument('-e','--entropy', default=0.0,type=float,help='upper cut-off for entropy, defaults to 0.0')
parser.add_argument('-u','--unknown', action='store_true',help='set -u, if contigs should be reported that were not annotated by Kraken')
parser.add_argument('-o','--outname', default="",help='name for the output can be specified or will be constructed from input file name and entropy cut-off')
parser.add_argument('-s','--silent', action='store_false',help='set -s, to suppress printing the number of uncontigs not annotated by Kraken')
parser.add_argument('-d','--taxdir', default=".",type=str,help='directory where the tax files sit')

args = parser.parse_args()
krakFile = args.inputFile
divthresh = args.entropy

taxdir = args.taxdir+"/"
if args.outname != "":
    outFile1 = args.outname
elif divthresh>0:
    outFile1 = krakFile + "annoEnt" + str(divthresh) + ".tsv"
else:
    outFile1 = krakFile + "annoUnambig.tsv"


# Definition of the class Node

class Node:
    """Noeud"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if tax_id in name_object:
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list

# Function to find common ancestor between two nodes or more
def common_ancestor(node_list):
    global name_object
    list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
    for node in node_list:
        list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
        ancestral_list = []                             
        for i in list1:
            if i in list2:                         # Identify common nodes between the two genealogy
                ancestral_list.append(i)                 
        list1 = ancestral_list                     # Reassing ancestral_list to list 1.
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestra_list is the common ancestor of all nodes.
    return common_ancestor                         # Return a node


#############################
#                           #
#   Read taxonomy files     #
#                           #
#############################

######################
# 
# Load names defintion

name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

# Load names file
namepath = taxdir + "all_tax2ID.tsv"
name_file =  open(namepath,"r")
while 1:
    line = name_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, name = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_reverse[name] = tax_id
    name_dict[tax_id] = name          # ... and load them into dictionary
name_file.close()

######################
# 
# Load kingdom definition

name_dict_kingdom = {}          # Initialise dictionary with TAX_ID:NAME

# Load kingdom file 
kingdompath = taxdir + "K_tax2ID.tsv"
kingdom_file =  open(kingdompath,"r")
while 1:
    line = kingdom_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, kingdom = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_kingdom[tax_id] = kingdom          # ... and load them into dictionary
kingdom_file.close()

######################
# 
# Load phylum definition

name_dict_phylum = {}          # Initialise dictionary with TAX_ID:NAME

# Load phylum file 
phylumpath = taxdir + "P_tax2ID.tsv"
phylum_file =  open(phylumpath,"r")
while 1:
    line = phylum_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, phylum = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_phylum[tax_id] = phylum          # ... and load them into dictionary
phylum_file.close()

######################
# 
# Load class definition

name_dict_pclass = {}          # Initialise dictionary with TAX_ID:NAME

# Load class file
pclasspath = taxdir + "C_tax2ID.tsv"
pclass_file =  open(pclasspath,"r")
while 1:
    line = pclass_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, pclass = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_pclass[tax_id] = pclass          # ... and load them into dictionary
pclass_file.close()


######################
# 
# Load order definition

name_dict_order = {}          # Initialise dictionary with TAX_ID:NAME

# Load order file made in R
orderpath = taxdir + "O_tax2ID.tsv"
order_file =  open(orderpath,"r")
while 1:
    line = order_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, order = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_order[tax_id] = order          # ... and load them into dictionary
order_file.close()

######################
# 
# Load family definition

name_dict_family = {}          # Initialise dictionary with TAX_ID:NAME

# Load family file 
familypath = taxdir + "F_tax2ID.tsv"
family_file =  open(familypath,"r")
while 1:
    line = family_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, family = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_family[tax_id] = family         # ... and load them into dictionary
family_file.close()

######################
# 
# Load genus definition

name_dict_genus = {}          # Initialise dictionary with TAX_ID:NAME

# Load genus file
genuspath = taxdir + "G_tax2ID.tsv"
genus_file =  open(genuspath,"r")
while 1:
    line = genus_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, genus = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_genus[tax_id] = genus         # ... and load them into dictionary
genus_file.close()

######################
# 
# Load species definition

name_dict_species = {}          # Initialise dictionary with TAX_ID:NAME

# Load species file 
speciespath = taxdir + "S_tax2ID.tsv"
species_file =  open(speciespath,"r")
while 1:
    line = species_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    tax_id, species = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_species[tax_id] = species         # ... and load them into dictionary
species_file.close()

######################
# 
# Load taxonomy

# Define taxonomy variable
global name_object
name_object = {}


# Load taxonomy NCBI file
taxonomypath = taxdir + "parents_IDs.tsv"
taxonomy_file = open(taxonomypath,"r")
while 1:
    line = taxonomy_file.readline()
    if line == "":
        break
    #print line
    line = line.rstrip()
    tab = line.split("\t")
    
    tax_id = str(tab[0])
    tax_id_parent = str(tab[1])

    # Define name of the taxid
    name = "unknown"
    if tax_id in name_dict:
        name = name_dict[tax_id]
    
    if tax_id not in name_object:
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id   = tax_id        # Assign tax_id
    name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
    name_object[tax_id].name     = name          # Assign name
    
    if  tax_id_parent in name_object:
        children = name_object[tax_id].children  # If parent is is already in the object
        children.append(tax_id)                  # ...we found its children.
        name_object[tax_id].children = children  # ... so add them to the parent
taxonomy_file.close()

#####################
#                   #
#contig annotation  #
#                   #
#####################
#function to calculate diversity of annotations
def shannonDiv(dictionary,sumTax):
    taxDiv = 0.0
    if len(dictionary) > 0 and sumTax > 1:
        for item in dictionary:
            taxDiv += (float(dictionary[item])/sumTax) * math.log(float(dictionary[item])/sumTax,2) / math.log(sumTax,2)
#            taxDiv += (float(dictionary[item])/sumTax) * math.log( float(dictionary[item])/sumTax)
#        taxDiv = taxDiv/math.log(sumTax)
    else:
        taxDiv = 1.0
    return 0.00 - taxDiv

# function to retrieve name and number of annotated bases for a taxon
def orgGenealCount(anc,taxDict,orgCnt,geneaDict,taxSum):
    if anc in geneaDict:
        taxName = geneaDict[anc]
        taxSum += int(orgCnt)
        if taxName not in taxDict:
            taxDict[taxName] = int(orgCnt)
        else:
            taxDict[taxName] += int(orgCnt)
    return taxDict, taxSum

# function to retrieve taxon name
def orgGenealName(anc,geneaDict,taxName):
    if anc in geneaDict:
        taxName = geneaDict[anc]
    return taxName

# function to test if lower taxon is in higher taxon and stop the annotation if necessary
def phyloTester(taxName,testList,retVal,annotationList):
    if taxName != "unknown":
        testList.append(name_dict_reverse[taxName])
        if len(testList) == 2:
            if testList[0] in name_object[testList[1]].genealogy():
                del testList[0]
            else:
                del annotationList[-1]
                retVal = -2.0
                taxName = "unknown"
                    

ucnt = 0

krak_file =  open(krakFile,"r")
out_file1 = open(outFile1, "w")
out_file1.write("contig" + "\t"+ "length" +"\t"+ "divSpec" + "\t"+"divGen" + "\t"+ "divFam" + "\t"+"divOrd" + "\t"+"divClass" + "\t"+"divPhylum" + "\t"+"divKingdom" + "\t"+ "annotationLevel" + "\t"+"species" + "\t" +"genus" + "\t"+"family" + "\t"+"order" + "\t"+"class" + "\t"+"phylum" + "\t"+"kingdom" +"\n")
while 1:
    linek = krak_file.readline()
    if linek == "":
        break
    linek = linek.rstrip()
    tabk = linek.split("\t")
    if tabk[0] == "U":
        if args.unknown:
            out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ "NA" + "\t"+ "NA" + "\t"+ "NA" + "\t"+ "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t"+ "not annotated" + "\t"+ "NA" + "\t" + "NA"  + "\t"+ "NA" + "\t"+"NA" + "\t"+"NA"+ "\t"+"NA" + "\t"+"NA" +"\n")
        ucnt += 1
    else:
        cotak = tabk[4].split(" ")
        orgs = {}
        specs = {}
        gens = {}
        fams = {}
        ords = {}
        clas = {}
        phys = {}
        kings = {}
        specsum = gensum = famsum = ordsum = clasum = physum = kingsum = 0
        for i in cotak:
            orgID = i.split(":")[0]
            orgCount = i.split(":")[1]
            if orgID not in orgs:
                orgs[orgID] = int(orgCount)
            else:
                orgs[orgID] += int(orgCount)
            if orgID in name_dict and orgID != "A" and orgID != "0":
                allOrg = []
                for anc in name_object[orgID].genealogy():
                    allOrg.append(anc)
                for anc in allOrg:
                    kings,kingsum = orgGenealCount(anc,kings,orgCount,name_dict_kingdom,kingsum)
                    phys,physum = orgGenealCount(anc,phys,orgCount,name_dict_phylum,physum)
                    clas,clasum = orgGenealCount(anc,clas,orgCount,name_dict_pclass,clasum)
                    ords,ordsum = orgGenealCount(anc,ords,orgCount,name_dict_order,ordsum)
                    fams,famsum = orgGenealCount(anc,fams,orgCount,name_dict_family,famsum)
                    gens,gensum = orgGenealCount(anc,gens,orgCount,name_dict_genus,gensum)
                    specs,specsum = orgGenealCount(anc,specs,orgCount,name_dict_species,specsum)
        specdiv = gendiv = famdiv = orddiv = cladiv = phydiv = kingdiv = -1.0
        specdiv = shannonDiv(specs,specsum)
        gendiv = shannonDiv(gens,gensum)
        famdiv = shannonDiv(fams,famsum)
        orddiv = shannonDiv(ords,ordsum)
        cladiv = shannonDiv(clas,clasum)
        phydiv = shannonDiv(phys,physum)
        kingdiv = shannonDiv(kings,kingsum)

        kingdom_nameB = phylum_nameB = pclass_nameB = order_nameB = family_nameB = genus_nameB = species_nameB = "unknown"
        if len(kings) >0:
            kingdom_nameB = max(kings.items(), key=operator.itemgetter(1))[0]
        if len(phys) >0:
            phylum_nameB = max(phys.items(), key=operator.itemgetter(1))[0]
        if len(clas) >0:
            pclass_nameB = max(clas.items(), key=operator.itemgetter(1))[0]
        if len(ords) >0:
            order_nameB = max(ords.items(), key=operator.itemgetter(1))[0]
        if len(fams) >0:
            family_nameB = max(fams.items(), key=operator.itemgetter(1))[0]
        if len(gens) >0:
            genus_nameB = max(gens.items(), key=operator.itemgetter(1))[0]
        if len(specs) >0:
            species_nameB = max(specs.items(), key=operator.itemgetter(1))[0]

        annoLev = ["none"]
        kingdom_nameA = phylum_nameA = pclass_nameA = order_nameA = family_nameA = genus_nameA = species_nameA = "unknown"
        testPhylo = []
        divthr = divthresh
        if kingdiv <= divthr:
            annoLev.append("kingdom")
            kingdom_nameA = kingdom_nameB
            if kingdom_nameA != "unknown":
                testPhylo.append(name_dict_reverse[kingdom_nameA])
            if phydiv <= divthr:
                annoLev.append("phylum")
                phylum_nameA = phylum_nameB
                phyloTester(phylum_nameA,testPhylo,divthr,annoLev)
                if cladiv <= divthr:
                    annoLev.append("class")
                    pclass_nameA = pclass_nameB
                    phyloTester(pclass_nameA,testPhylo,divthr,annoLev)
                    if orddiv <= divthr:
                        annoLev.append("order")
                        order_nameA = order_nameB
                        phyloTester(order_nameA,testPhylo,divthr,annoLev)
                        if famdiv <= divthr:
                            annoLev.append("family")
                            family_nameA = family_nameB
                            phyloTester(family_nameA,testPhylo,divthr,annoLev)
                            if gendiv <= divthr:
                                annoLev.append("genus")
                                genus_nameA = genus_nameB
                                phyloTester(genus_nameA,testPhylo,divthr,annoLev)
                                if specdiv <= divthr:
                                    annoLev.append("species")
                                    species_nameA = species_nameB
                                    phyloTester(species_nameA,testPhylo,divthr,annoLev)
        annoLevel = annoLev[-1]
        if annoLevel == "species" and species_nameA == "unknown":
            annoLevel = "genus"
        if annoLevel == "genus" and genus_nameA == "unknown":
            annoLevel = "family"
        if annoLevel == "family" and family_nameA == "unknown":
            annoLevel = "order"
        if annoLevel == "order" and order_nameA == "unknown":
            annoLevel = "class"
        if annoLevel == "class" and pclass_nameA == "unknown":
            annoLevel = "phylum"
        if annoLevel == "phylum" and phylum_nameA == "unknown":
            annoLevel = "kingdom"
        if annoLevel == "kingdom" and kingdom_nameA == "unknown":
            annoLevel = "none"
        
        out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ str(specdiv) + "\t"+ str(gendiv) + "\t"+ str(famdiv) + "\t"+ str(orddiv) + "\t" + str(cladiv) + "\t" + str(phydiv) + "\t" + str(kingdiv) + "\t"+ annoLevel + "\t"+species_nameA + "\t" +genus_nameA + "\t"+family_nameA + "\t"+order_nameA + "\t"+pclass_nameA+ "\t"+phylum_nameA + "\t"+kingdom_nameA +"\n")
        
krak_file.close()
out_file1.close()

if args.silent:
    print(ucnt) 

