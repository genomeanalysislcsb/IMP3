#!/usr/bin/env python3

# Script to filter the gff annotation file from the filtered .fa assembly file 
# Get the contigs of the annotation file that exist in the assembly

#import library
import sys

# read input files
ann_file = snakemake.input.in_gff
ass_file = snakemake.input.in_ass
    
#Open the files in reading mode
annotations = open(ann_file, "r")
assemblies = open(ass_file, "r")

#open output file in writing mode
new_ann_file = snakemake.output[0]
new_ann_handle = open(new_ann_file,'w')
#Results in a list of lines, in column organization 
new_annotations = [] 


contig_IDs = []
#Make a list of the contig names in the assembly file by iterating over the lines in the assembly file
for line in assemblies:
    #check if the line starts with > symbol, indicating the location of a contig name
    if line.startswith('>'):
        #assign to variable
        contig_ID = line 
        #remove the '>' and the '\n' from the contigID 
        contig_ID = contig_ID.replace('>','')
        contig_ID = contig_ID.replace('\n','')
        #add to a list of contigs IDs
        contig_IDs.append(contig_ID) 

#iterate over lines in annotation file 
for line in annotations:
    #Separate the columns from the .gff file and take first element as value variable
    columns = line.strip().split('\t')
    #obtain the first column as it obtains the contigID
    value = columns[0]
    #check if that contigID is in the list from the assembly file
    if value in contig_IDs:
        #add the line to the file 
        new_annotations.append(line)
    else:
        #pass all the contigs which are not in the assembly file, basically filtering them out
        pass 

#Write all the found contigs into the output file 
for element in new_annotations:
    new_ann_handle.write(element + '\n')


#Remember to close files after done
annotations.close()
assemblies.close() 
new_ann_handle.close()
