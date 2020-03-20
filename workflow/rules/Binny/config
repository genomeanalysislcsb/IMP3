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

def dict_merge(a, b):
    """
    Deep merge 2 dicts together
    """
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
            result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result


# default configuration file
configfile:
    srcdir("../../../config/config.binny.yaml")


# default executable for snakmake
shell.executable("bash")

# custom configuration file
CUSTOM_CONFIG_PATH = os.environ.get("CONFIGFILE", "../../../config/userconf.binny.yaml")

# merge 2 configurations files together
if os.path.exists(CUSTOM_CONFIG_PATH):
    print("reading configuration from "+CUSTOM_CONFIG_PATH)
    with open(CUSTOM_CONFIG_PATH, 'r') as rhandle:
        data = yaml.load(rhandle)
        config = dict_merge(config, data)


# some parameters
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")

# get parameters from the command line
# output
OUTPUTDIR = os.environ.get("OUTPUTDIR", config['outputdir'])
SUMDIR = os.environ.get("SUMDIR", config['summarydir'])
if OUTPUTDIR == SUMDIR:
    if SUMDIR[-1] != "/":
        SUMDIR += "_final"
    else:
        SUMDIR = SUMDIR[:-1] + "_final"

# input
CONTIGS = os.environ.get("CONTIGS", config['raws']['Contigs'])
MGaln = os.environ.get("MGaln", config['raws']['Alignment_metagenomics'])
SAMPLE = os.environ.get("SAMPLE", config['sample'])
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.environ.get("DBPATH", os.path.expandvars(config['db_path']))
if not os.path.exists(DBPATH):
    os.makedirs(DBPATH)
    
# define the data types used and the assembly
TYPES = ['mg']
ASS = 'mg'

# settings
BINNY_STEPS = ["analysis","binning"]  
BIN_STEPS = ["binny"]
SUM_STEPS = os.environ.get("SUM_STEPS", config['summary_steps']).split()

# hardware parameters
BIGCORENO = os.environ.get("BIGCORENO", config['mem']['big_mem_cores'])
BIGMEMTOTAL = os.environ.get("BIGMEMTOTAL", config['mem']['big_mem_total_gb'])
MEMCORE = os.environ.get("MEMCORE", str(config['mem']['normal_mem_per_core_gb']) + "G")
BIGMEMCORE = os.environ.get("BIGMEMCORE", str(config['mem']['big_mem_per_core_gb']) + "G")

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = os.environ.get("TMPDIR", config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

