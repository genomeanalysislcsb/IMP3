#!/usr/bin/env python

import os
import shutil
import subprocess

if len(snakemake.input) != len(snakemake.output):
    raise OSError("//Inputs and outputs are not of the same length: %s <> %s" % (', '.join(snakemake.input), ', '.join(snakemake.output)))

for ifile, ofile in zip(snakemake.input, snakemake.output):
    print(ifile, '=>', ofile)
    # ungunzip
    if os.path.splitext(ifile)[-1] in ['.gz', '.gzip']:
        subprocess.check_call("pigz -dc -p {p} {i} > {o}".format(p=snakemake.threads, i=ifile, o=ofile), shell=True)
    # unbzip2
    elif os.path.splitext(ifile)[-1] in ['.bz2', '.bzip2']:
        subprocess.check_call("pbzip2 -dc -p {p} {i} > {o}".format(p=snakemake.threads, i=ifile, o=ofile), shell=True)
    # copy
    else:
        shutil.copy(ifile, ofile)