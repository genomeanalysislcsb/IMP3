#!/usr/bin/env python3

import os
import sys
import argparse
import operator
import math
import numpy as np

parser = argparse.ArgumentParser(description='Process output of Kraken run on contigs.')
parser.add_argument('-c','--contig', default="", help='per-contig output from Kraken')
parser.add_argument('-r','--report', default="", help='report output from Kraken (all taxa must be reported)')
parser.add_argument('-e','--entropy', default=0.0,type=float,help='upper cut-off for entropy, defaults to 0.0')
parser.add_argument('-u','--unknown', action='store_true',help='set -u, if contigs should be reported that were not annotated by Kraken')
parser.add_argument('-o','--outname', default="",help='name for the output can be specified or will be constructed from input file name and entropy cut-off')
parser.add_argument('-s','--silent', action='store_false',help='set -s, to suppress printing the number of contigs not annotated by Kraken')
#parser.add_argument('-d','--taxdir', default=".",type=str,help='directory where the tax files sit')

args = parser.parse_args()
krakFile = args.contig
taxFile = args.report
divthresh = args.entropy

#taxdir = args.taxdir+"/"
if args.outname != "":
    outFile1 = args.outname
elif divthresh>0:
    outFile1 = krakFile + "annoEnt" + str(divthresh) + ".tsv"
else:
    outFile1 = krakFile + "annoUnambig.tsv"


# Definition of the class Node

class Node:
    """Node"""
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
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestral_list is the common ancestor of all nodes.
    return common_ancestor                         # Return a node


#############################
#                           #
#   Read taxonomy file      #
#                           #
#############################

global name_object
name_object = {}
name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID #not all names are unique, though
parentList = [None] * 100
rank_dict = {}          # Initialise dictionary with TAX_ID:RANK
rank_dict_reverse = {}  # Initialise dictionary with RANK:[TAX_IDs]
rank_par_dict = {}      # dictionary with RANK:parent_RANK
rank_child_dict = {}    # dictionary with RANK:child_RANK

tax_file = open(taxFile, "r")
old_ilevel = 0
while 1:
    line = tax_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")   # 0: total perc, 1:total counts, 2: counts at this level, 3: rank abbrev., 4: tax ID, 5: name (indented)
    if tab[3] == "U":
        next
    rank, tax_id, indentName = tab[3], tab[4], tab[5]      # Assign tax_id and name ...
    name = indentName.lstrip(' ')
    ilevel = int((len(indentName) - len(name))/2)     # get level from indentation
    print(ilevel)
    parentList[ilevel] = tax_id
    if ilevel > 0:
        tax_id_parent = parentList[ilevel-1]
    if name not in name_dict_reverse:
        name_dict_reverse[name] = tax_id
        name_dict[tax_id] = name                # ... and load them into dictionary
    else:
        if ilevel > 0:
            name_dict_reverse[name + "." + name_dict[tax_id_parent]] = tax_id
            name_dict[tax_id] = name  + "." + name_dict[tax_id_parent]  # ... for kids that have their parent's names or worse
                        
    rank_dict[tax_id] = rank
    if rank not in rank_dict_reverse:
        rank_dict_reverse[rank] = [tax_id]
    else:
        rank_dict_reverse[rank].append(tax_id)
        
    if ilevel > 0:                                                  #fill dictionaries of rank IDs to reconstruct taxonomy tree later
        if rank not in rank_par_dict:
            rank_par_dict[rank] = {rank_dict[tax_id_parent]: 1}
        elif rank_dict[tax_id_parent] not in rank_par_dict[rank]:
            rank_par_dict[rank][rank_dict[tax_id_parent]] = 1
        if rank_dict[tax_id_parent] not in rank_child_dict:
            rank_child_dict[rank_dict[tax_id_parent]] = {rank: 1}
        elif rank not in rank_child_dict[rank_dict[tax_id_parent]]:
            rank_child_dict[rank_dict[tax_id_parent]][rank] = 1
        
    if tax_id not in name_object:          # this should always be the case - otherwise we would overwrite in the following lines
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id = tax_id             # Assign tax_id
    if ilevel > 0:
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
    name_object[tax_id].name = name                 # Assign name
    
    if ilevel > 0:
        siblings = name_object[tax_id_parent].children  # Parent is always already in the object
        siblings.append(tax_id)                         # ...we found its children.
        name_object[tax_id_parent].children = siblings  # ... so add them to the parent 
tax_file.close()

########################
#                      #
# reconstruct taxonomy #
#                      #
########################

ranklist = [rank_dict["1"]]
while 1:
 #   print(" ".join(ranklist))
    if ranklist[-1] not in rank_child_dict:
        break
    pot_child = list(rank_child_dict[ranklist[-1]].keys())   # find potential next ranks in ordered list
 #   print(" ".join(pot_child))
    if len(pot_child) == 1:                           # if there's only one child, that's it
        ranklist.append(pot_child[0])
    elif len(pot_child) > 1:                          # if there are more than one
        keep_child = []
        for cchild in pot_child:                      # check all potential next ranks
 #           print(" ".join(rank_par_dict[cchild].keys()))
            if set(rank_par_dict[cchild].keys()).isdisjoint(set(pot_child)):      # if they're children of other candidates, remove
                keep_child.append(cchild)
 #       print(" ".join(keep_child))
        if len(keep_child) == 1:                      # remaining candidates can be 1 or more
            ranklist.append(keep_child[0])            # if one, this is the next one
        else:                                         
            keep_par = []                             #otherwise, check which of the candidates is a parent of other candidates
            for pot_par in keep_child:
                if not set(rank_child_dict[pot_par].keys()).isdisjoint(set(pot_child)):
                    keep_par.append(pot_par)
            if len(keep_par) == 1:                    # if there's just one parent , this one is next
                ranklist.append(keep_par[0])
            elif len(keep_par) > 1:                   # if more than one parent of other candidates, check if one of them is a child of an earlier rank
                keep_grandchild = []
                for poty in keep_par:
                    if not set(ranklist[:-1]).isdisjoint(set(rank_par_dict[poty].keys())):
                        keep_grandchild.append(poty)
                if len(keep_grandchild) == 1:
                    ranklist.append(keep_grandchild[0])     # take the child of an earlier rank (at least it's not an in-between rank)
                else:                                  #we have no more options to establish order based on relationships
                    keep_letter = []                   # we try based on the letter
                    start_par = ranklist[-1][0]
                    for potx in keep_par:
                        start_pot = potx[0]
                        if start_pot != start_par:
                            keep_letter.append(potx)   # keeping the one that isn't the same letter as the parent
                    if len(keep_letter) == 1:            
                        ranklist.append(keep_letter[0])
                    else:                                # if there's more than one candidate, we're dead
                        sys.exit("Taxonomy cannot be inferred based on the data.")
            else:                                  #we have no more options to establish order based on relationships
                print(" ".join(keep_child))
                keep_letter = []                   # we try based on the letter
                start_par = ranklist[-1][0]
                for potx in keep_child:
                    start_pot = potx[0]
                    if start_pot != start_par:
                        keep_letter.append(potx)   # keeping the one that isn't the same letter as the parent
                if len(keep_letter) == 1:            
                    ranklist.append(keep_letter[0])
                else:                                # if there's more than one candidate, we're dead
                    sys.exit("Taxonomy cannot be inferred based on the data.")
                        
    else:                                            # if there are no more children, we've reached the end             
        break                                        # we still need to add more branches for the in-between ranks
i = 0
while i < len(ranklist):
    if not ranklist[i] in rank_child_dict:
        break
    missing_son = set(rank_child_dict[ranklist[i]].keys()) - set(ranklist)
 #   print(missing_son)
    if len(missing_son) > 0:      #if one child has not been added yet, we find it
        if len(missing_son) != 1:
            sys.exit("Something is wrong.")
        inserter = missing_son.pop()
        ranklist.insert(i+1,inserter)            # and insert it right after the parent
    i += 1

print(" ".join(ranklist))
   

#####################
#                   #
# contig annotation #
#                   #
#####################
#function to calculate diversity of annotations
def shannonDiv(dictionary,sumTax):
    taxDiv = 0.0
    if len(dictionary) > 0 and sumTax > 1:
        for item in dictionary:
            taxDiv += (float(dictionary[item])/sumTax) * math.log(float(dictionary[item])/sumTax,2) / math.log(sumTax,2)
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
out_file1.write("contig" + "\t"+ "length" +"\t"+ "entropy" + "\t"+ "annotationLevel" + "\t"+ "\t".join(ranklist) +"\n")
while 1:
    linek = krak_file.readline()
    if linek == "":
        break
    linek = linek.rstrip()
    tabk = linek.split("\t")
    if tabk[0] == "U":
        if args.unknown:
            out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ "NA" + "\t"+ "not annotated" + "\t"+ "\t".join(["NA"] * len(ranklist)) +"\n")
        ucnt += 1
    else:
        cotak = tabk[4].split(" ")
        orgs = {}
        for i in cotak:
            orgID = i.split(":")[0]
            orgCount = i.split(":")[1]
            if orgID not in orgs:
                orgs[orgID] = int(orgCount)
            else:
                orgs[orgID] += int(orgCount)
            if orgID in name_dict and orgID != "A" and orgID != "0":
                for anc in name_object[orgID].genealogy(): #should return taxids
                    if anc not in orgs:
                        orgs[anc] = int(orgCount)
                    else:
                        orgs[anc] += int(orgCount)
        ranksum_dict = {}
        rank_orgs = {}
        for ctax in orgs.keys():
            if rank_dict[ctax] not in ranksum_dict:
                ranksum_dict[rank_dict[ctax]] = orgs[ctax]
                rank_orgs[rank_dict[ctax]] = {ctax: orgs[ctax]}
            else:
                ranksum_dict[rank_dict[ctax]] += orgs[ctax]
                rank_orgs[rank_dict[ctax]][ctax] = orgs[ctax]
        cranks = sorted(ranksum_dict.items(), key=operator.itemgetter(1), reverse = True)
        rankwin_dict = {}
        winnerList = []
        first = True
        for rank in cranks:
            cdiv = shannonDiv(rank_orgs[rank],ranksum_dict[rank])
            if cdiv <= divthr:
                cwin = max(rank_orgs[rank].items(), key=operator.itemgetter(1))[0]
                winancs = name_object[cwin].genealogy()
                if first:
                    rankwin_dict[rank] = cwin
                    winnerList.append(cwin)
                    first = False
                elif not set(winnerList).isdisjoint(set(winancs)):
                    rankwin_dict[rank] = cwin
                    winnerList.append(cwin)
            else:
                break
        #out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ "NA" + "\t"+ "not annotated" + "\t"+ "\t".join(["NA"] * len(ranklist)) +"\n")
        ranked_winner_list = []
        went = "NA"
        wrank = "not annotated"
        for crank in ranklist:
            if crank in rankwin_dict:
                ranked_winner_list.append(name_dict[rankwin_dict[crank]])
                went = shannonDiv(rank_orgs[crank],ranksum_dict[crank])
                wrank = crank
            else:
                ranked_winner_list.append("NA")
        if went != "NA" or args.unknown:
            out_file1.write(tabk[1] + "\t"+ tabk[3] +"\t"+ went + "\t"+ wrank + "\t"+ "\t".join(ranked_winner_list) +"\n")        
krak_file.close()
out_file1.close()

if args.silent:
    print(ucnt)

