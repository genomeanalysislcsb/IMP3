import os
import sys
import shutil
import gzip
#import json
import yaml
import bz2
import re
from copy import deepcopy
import subprocess
import pandas as pd
from pathlib import Path
from snakemake.utils import validate

# NOTE: not used
# def dict_merge(a, b):
#     """
#     Deep merge 2 dicts together
#     """
#     if not isinstance(b, dict):
#         return b
#     result = deepcopy(a)
#     for k, v in b.items():
#         if k in result and isinstance(result[k], dict):
#             result[k] = dict_merge(result[k], v)
#         else:
#             result[k] = deepcopy(v)
#     return result

# default configuration file
configfile:
    srcdir("../../../config/config.imp.yaml")


# default executable for snakmake
shell.executable("bash")

# NOTE: not used (rm os.environ.get() and path is not valid)
# # custom configuration file
# CUSTOM_CONFIG_PATH = os.environ.get("CONFIGFILE", "../../../config/userconf.imp.yaml")
# # merge 2 configurations files together
# if os.path.exists(CUSTOM_CONFIG_PATH):
#     #print("reading configuration from "+CUSTOM_CONFIG_PATH)
#     with open(CUSTOM_CONFIG_PATH, 'r') as rhandle:
#         data = yaml.load(rhandle)
#         config = dict_merge(config, data)

# config validation
validate(config,schema="../../../schemas/schema.imp.yaml")

# relevant paths
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")
MODDIR = srcdir("../../../submodules") # git submodules

# output
OUTPUTDIR = os.path.abspath(os.path.expandvars(config['outputdir']))
SUMDIR    = os.path.abspath(os.path.expandvars(config['summarydir']))
if OUTPUTDIR == SUMDIR:
    if SUMDIR[-1] != "/":
        SUMDIR += "_final"
    else:
        SUMDIR = SUMDIR[:-1] + "_final"

# input
MG      = config['raws']['Metagenomics'].split()
MT      = config['raws']['Metatranscriptomics'].split()
LONG    = config['raws']['LongReads']
CONTIGS = config['raws']['Contigs']
MGaln   = config['raws']['Alignment_metagenomics']
MTaln   = config['raws']['Alignment_metatranscriptomics']
ANN     = config['raws']['Gff']
SAMPLE  = config['sample']
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.path.expandvars(config['db_path'])
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)
    
# define the data types used and the assembly
if MG and MT and config['assembly']['hybrid']:
    TYPES = ['mg', 'mt']
    ASS = 'mgmt'
elif MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif MG:
    TYPES = ['mg']
    ASS = 'mg'
elif MT:
    TYPES = ['mt']
    ASS = 'mt'
elif MGaln and MTaln:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif MGaln:
    TYPES = ['mg']
    ASS = 'mg'
elif MTaln:
    TYPES = ['mt']
    ASS = 'mt'
elif CONTIGS and MG and MT:
    TYPES = ['mg', 'mt']
    ASS = 'mg'
elif CONTIGS and MG:
    TYPES = ['mg']
    ASS = 'mg'
elif CONTIGS and MT:
    TYPES = ['mt']
    ASS = 'mg'

#print(TYPES)

# steps
IMP_STEPS               = sorted(config['steps'].split())
PREPROCESSING_FILTERING = config['preprocessing_filtering']
FILTER                  = config['filtering']['filter'].split()
IMP_ASSEMBLER           = config['assembly']['assembler']
SUM_STEPS               = config['summary_steps'].split()
BIN_STEPS               = config['binning']['binners'].split()
if not MG and not MGaln:
    BIN_STEPS = list(set(BIN_STEPS).difference(["MaxBin", "binny"]))

# hardware parameters
BIGCORENO   = config['mem']['big_mem_cores']
BIGMEMTOTAL = config['mem']['big_mem_total_gb']
MEMCORE     = str(config['mem']['normal_mem_per_core_gb']) + "G"
BIGMEMCORE  = str(config['mem']['big_mem_per_core_gb']) + "G"

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = config['tmp_dir']
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)
