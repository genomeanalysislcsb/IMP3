#!/usr/bin/env python3

#Script to parse all regions Genbank files from antismash output to a single gff file 
#also add the MiBig identifiers that have the hoghest scoring hits from the knownclusterblastoutput.txt file
#and add it to the additional attributes 

#import libraries 
import glob 
import os
import sys 
from Bio import SeqIO


#take the input argument - directory path to where all the files to read are 
input_path = sys.argv[1]
#The glob module finds all the pathnames matching a specified pattern
#list of files from directory folder 
gbk_glob = os.path.join(input_path, "*contig*.gbk")
filenames = glob.glob(gbk_glob)
cb_txtfile = os.path.join(input_path, "knownclusterblastoutput.txt") 



#function to write the genbank file region features into lines for each file(contig)
def make_lines(file_path):

    #read each file as a SeqIO object with file type as genbank
    genome_record = SeqIO.read(file_path, 'genbank')
    

    #make a list to store each line as an element
    lines = []
    line = ''
    #iterate over all the features(sections of the contig) annotated by the tool 
    for i in range(len(genome_record.features)):


        #different types of features have different components(qualifiers/keys) - check to avoid errors when writting attributes for each type
        #check if the type of region/feature is a coding sequence/annotated region with the 'phase' qualifier/key, to also include along the line
        if genome_record.features[i].qualifiers.__contains__('phase'):
            #check if annotation is from a 'tool' a 'aStool' or a 'source' to include in the second column 
            if genome_record.features[i].qualifiers.__contains__('tool'):
                string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['tool'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(genome_record.features[i].qualifiers['phase'][0]) + '\t'
            else:
                if genome_record.features[i].qualifiers.__contains__('aStool'):
                    string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['aStool'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(genome_record.features[i].qualifiers['phase'][0]) + '\t'
                else:
                    string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['source'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(genome_record.features[i].qualifiers['phase'][0]) + '\t'
        else:
            if genome_record.features[i].qualifiers.__contains__('tool'):
                string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['tool'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(0) + '\t' 
            else:
                if genome_record.features[i].qualifiers.__contains__('aStool'):
                    string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['aStool'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(0) + '\t'
                else:
                    string = str(genome_record.description) + '\t' + str(genome_record.features[i].qualifiers['source'][0]) + '\t' + str(genome_record.features[i].type) + '\t' + str(genome_record.features[i].location.start + 1) + '\t' + str(genome_record.features[i].location.end) + '\t' + '.' + '\t' + str(genome_record.features[i].location)[-2] + '\t' + str(0) + '\t'
        

        #collect all additional information of each line as one string
        attributess = ''
        #iterate over the names of the qualifiers, which are keys in a dictionary 
        for key in genome_record.features[i].qualifiers:
            #exclude the qualifiers with the following key names, as they are redundant 
            if key ==  'ID' or key == 'gene' or key == 'translation' or key == 'partial' or key == 'domain_id':
                pass
            else:
                #check if a certain qualifier has a values list larger than one, then include all 
                if len(genome_record.features[i].qualifiers[key]) > 1:
                    temp = str(key) + '=' + str(genome_record.features[i].qualifiers[key]) + ';' 
                    #save the attributes in the string one by one as they iterate
                    attributess = attributess + str(temp)
                else:
                    #if not, select only the one value, which is the first element of the list
                    temp = str(key) + '=' + str(genome_record.features[i].qualifiers[key][0]) + ';' 
                    #save the attributes in the string one by one as they iterate
                    attributess = attributess + str(temp)

        #add all the component to make the complete line for each annotation 
        line = line + string + attributess
        if line.__contains__(';;'):
            line = line.replace(';;',';')
        #add each line as an element of the list 
        lines.append(str(line))
        #clear line for next iteration
        line = ''
    #return list of all the lines 
    return lines



#loop over all the gbk files in the list and apply the above function to all
file_list = []
#for each file(by calling the name) in the list of file names
for file in filenames:
    tmp_file = make_lines(file)
    file_list = file_list + tmp_file



#open knownclusterblast file from list 
f = open(cb_txtfile, 'r')


#dictionary to store the MiBig identifier with the corresponding locus_tag 
locus_dict = {}    

#iterate over the lines in the knownclusterblast text file
for Line in f:
    #lines are separated by | symbol, so split into columns
    columns = Line.split('|')
    #column with index 4 is the locus_tag 
    locus = columns[4]
    #column 5 contains 2 information elements separated by a tab, split and grab the second element as that is the MiBig identifier
    BGC_id = columns[5].split('\t')[1]
    #The highest is the first so if the key is already in the dictionary skip so no others are added for that locus tag
    if locus in locus_dict.keys():
        pass
    else:
        locus_dict[locus] = BGC_id

#iterate over the dictionary, separating both item types 
for key, value in locus_dict.items():
    #iterate over the lines list, enumerate the index as well as the actual content of the line
    for index, line in enumerate(file_list):
        #check if the dictionary key - locus_tag - is in the text file line 
        if key in line:
            #add the MiBig identifier to the line
            line = line + 'MiBig=' + value + ';'
            #update the lines list 
            file_list[index] = line
             
#close text file 
f.close()

    

#write gff outfile with all the types entrees for antismash annotation from gbk file 
all_out = os.path.join(input_path, "antismash.all.gff")
out_f = open(all_out, 'w')
out_f.write('\n'.join(file_list))
out_f.close()

#output files
gf_out = os.path.join(input_path, "antismash.gene_functions.gff")
gf_out_f = open(gf_out,'w')
pfam_out = os.path.join(input_path, "antismash.db_xref.gff")
pfam_out_f = open(pfam_out,'w')
aSDomain_out = os.path.join(input_path, "antismash.aSDomain.gff")
aSDomain_out_f = open(aSDomain_out, 'w')
mibig_out = os.path.join(input_path, "antismash.MiBig.gff")
mibig_out_f = open(mibig_out, 'w')

for line in file_list:
#write gff outfile for CDS regions which have actual annotated locus tags with identified gene_functions
    if line.__contains__('gene_functions'):
        gf_out_f.write(line + '\n')
#write gff outfile for found regions from Pfam database 
    if line.__contains__('PFAM_domain'):
        line = line.replace('PFAM_domain', 'CDS')
        pfam_out_f.write(line + '\n')
#write gff outfile for aSDomain regions 
    if line.__contains__('	aSDomain	'):
        line = line.replace('	aSDomain	', '	CDS	')
        aSDomain_out_f.write(line + '\n')
#write gff outfile for regions with MiBig annotations 
    if line.__contains__('MiBig'):
        if line.__contains__('; E-value:'):
            line = line.replace('; E-value:',': E-value')
        elif line.__contains__('PFAM_domain'):
            line = line.replace('PFAM_domain','CDS')
        elif line.__contains__('	aSDomain	'):
            line = line.replace('	aSDomain	', '	CDS	')
        elif line.__contains__('aSModule'):
            line = line.replace('aSModule', 'CDS')
        elif line.__contains__('CDS_motif'):
            line = line.replace('CDS_motif', 'CDS')
        mibig_out_f.write(line + '\n')

#close all files when done with them
gf_out_f.close()
pfam_out_f.close() 
aSDomain_out_f.close()
mibig_out_f.close()
