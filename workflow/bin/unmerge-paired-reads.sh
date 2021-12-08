#!/bin/bash
#
#	script: unmerges interweaving paired-reads into two separate files 
#
# (1) Scenario: Paired-end forward and reverse reads interweaved by merge-paired-reads.sh
#
#          READ 1
#    |-------------->
# 5' |-------------------------------------------------| 3'
#      | | | | | | | | | | | | | | | | | | | | | | | |
# 3' |-------------------------------------------------| 5'
#                                       <--------------|
#                                            READ 2
#
#    See "Example 5: sortmerna on forward-reverse paired-end
#    reads (2 input files)" of the SortMeRNA User Manual (version 1.7 
#    and higher)
#
#    Use unmerge-paired-reads.sh to separate the interweaved reads from
#    a single file into two separate files, where READ 1 will be
#    on the same line as READ 2 but in two separate files
#
# command: bash unmerge-paired-reads.sh inputfile.fastq file1.fastq file2.fastq 
#
#
# date: May 13, 2013
# contact: evguenia.kopylova@lifl.fr
#
#
# GPL-3.0 License 
# script is part of sortmeRNA: https://github.com/biocore/sortmerna
# please cite:
# Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611
#    
#    (c) evguenia kopylova 2013
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



# check all files are given
if [ $# -lt 3 ]; then
 echo "usage: $0 merged-reads-file.fastq forward-reads-name.fastq reverse-reads-name.fastq"
 exit 2
elif [ $# -gt 3 ]; then
 echo "usage: $0 merged-reads-file.fastq forward-reads-name.fastq reverse-reads-name.fastq "
 exit 2
fi

# quickly check input file is not fasta
sig1="$(head -c 1 $1)"
if [ "$sig1" == ">" ]; then
 echo "   Warning: $1 seems to be in fasta format (this script is for fastq reads only)"
fi


# unmerge paired reads
echo "   Processing $2 .."
perl -pe 's/\n/\t/ if $. %4' "$1" | awk 'NR%2 {print}' | tr "\t" "\n" >| "$2"
echo "   Processing $3 .."
perl -pe 's/\n/\t/ if $. %4' "$1" | awk '(NR+1)%2 {print}' | tr "\t" "\n" >| "$3"

echo "   Done."


