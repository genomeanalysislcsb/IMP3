#!/usr/bin/env python3
import sys
import textwrap

dic = {}
tmp_contig = ""
with open(sys.argv[1]) as file:  # Use with statement for better file handling
    for line in file:
        if line[0] == ">":
            if tmp_contig != "":
                dic[name] = tmp_contig
                tmp_contig = ""
            name = line.strip()
        else:
            tmp_contig += line.strip()
dic[name] = tmp_contig

for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
    print(k)  # Update print syntax for Python 3
    print(textwrap.fill(dic[k], 100, break_on_hyphens=False))
