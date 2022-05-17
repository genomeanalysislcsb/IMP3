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


# default configuration file
configfile:
    srcdir("../../../config/config.imp.yaml")


# default executable for snakmake
shell.executable("bash")

validate(config,schema="../../../schemas/schema.imp.yaml")


if config['sessionKind'] == "cluster":
    if str(config['settingsLocked']).lower() in ['true','t']:
        config['settingsLocked'] = True
        ini_normalMem = config['locked_normalMem']
        ini_bigMem = config['locked_bigMem']
        ini_bigCores = config['locked_bigCores']
        ini_bigTotal = config['locked_bigTotal']

        if ini_normalMem == 0 or ini_bigMem == 0 or ini_bigCores == 0:
            raise Exception("You're attempting to submit IMP3 to a cluster and but you've requested 0 cores or 0G memory per core. Ask the person who set up IMP3 to specify proper values in " + ROOTDIR + "/VARIABLE_CONFIG")
        config['mem']['normal_mem_per_core_gb'] = int(ini_normalMem)
        config['mem']['big_mem_per_core_gb'] = int(ini_bigMem)
        config['mem']['big_mem_cores'] = int(ini_bigCores)
        config['mem']['big_mem_total_gb'] = max(int(ini_bigTotal), int(config['mem']['big_mem_cores']) * int(config['mem']['big_mem_per_core_gb']))
        print("These memory requirements will be used:")
        print("normal mem: " + str(config['mem']['normal_mem_per_core_gb']) + "G")
        print("big mem: " + str(config['mem']['big_mem_per_core_gb']) + "G")
        print("big mem cores: " + str(config['mem']['big_mem_cores']))
        print("big mem total: " + str(config['mem']['big_mem_total_gb']) + "G")
    else:
        config['settingsLocked'] = False
elif config['sessionKind'] == "dryrun":
    if str(config['settingsLocked']).lower() in ['true','t']:
        config['settingsLocked'] = True
        ini_normalMem = config['locked_normalMem']
        ini_bigMem = config['locked_bigMem']
        ini_bigCores = config['locked_bigCores']
        ini_bigTotal = config['locked_bigTotal']
        config['mem']['normal_mem_per_core_gb'] = int(ini_normalMem)
        config['mem']['big_mem_per_core_gb'] = int(ini_bigMem)
        config['mem']['big_mem_cores'] = int(ini_bigCores)
        config['mem']['big_mem_total_gb'] = max(int(ini_bigTotal), int(config['mem']['big_mem_cores']) * int(config['mem']['big_mem_per_core_gb']))
        if ini_normalMem == 0 or ini_bigMem == 0 or ini_bigCores == 0:
            print("You've requested 0 cores or 0G memory per core. Ask the person who set up IMP3 to specify proper values in " + ROOTDIR + "/VARIABLE_CONFIG")
        else: 
            print("These memory requirements would be used:")
            print("normal mem: " + str(config['mem']['normal_mem_per_core_gb']) + "G")
            print("big mem: " + str(config['mem']['big_mem_per_core_gb']) + "G")
            print("big mem cores: " + str(config['mem']['big_mem_cores']))
            print("big mem total: " + str(config['mem']['big_mem_total_gb']) + "G")
    else:
        config['settingsLocked'] = False
else:
    config['mem']['big_mem_cores'] = workflow.cores
    if(config['mem']['big_mem_per_core_gb'] != config['mem']['normal_mem_per_core_gb']):
        config['mem']['normal_mem_per_core_gb'] = max(config['mem']['normal_mem_per_core_gb'],config['mem']['big_mem_per_core_gb'])
        config['mem']['big_mem_per_core_gb'] = max(config['mem']['normal_mem_per_core_gb'],config['mem']['big_mem_per_core_gb'])
        config['mem']['big_mem_total_gb'] = config['mem']['big_mem_cores'] * config['mem']['big_mem_per_core_gb']
        print("adjusted memory settings because there is only one kind of core in this run. CAUTION: we're using the bigger number of normal_mem_per_core_gb and big_mem_per_core_gb and the number of threads granted to IMP3.")
    print("Final resource settings:")
    print("maximal number of cores: " + str(workflow.cores))
    print("normal mem: " + str(config['mem']['normal_mem_per_core_gb']) + "G")
    print("big mem: " + str(config['mem']['big_mem_per_core_gb']) + "G")
    print("big mem cores: " + str(config['mem']['big_mem_cores']))
    print("big mem total: " + str(config['mem']['big_mem_total_gb']) + "G")


# some parameters
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")

# get parameters from the command line
# output
OUTPUTDIR = os.path.abspath(os.path.expandvars(config['outputdir']))
#print(OUTPUTDIR)
SUMDIR = os.path.abspath(os.path.expandvars(config['summarydir']))
if OUTPUTDIR == SUMDIR:
    if SUMDIR[-1] != "/":
        SUMDIR += "_final"
    else:
        SUMDIR = SUMDIR[:-1] + "_final"

EMAIL = config['email']
if EMAIL != "":
    if not re.fullmatch(r"^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$", EMAIL):
        EMAIL = ""
        print("Your email address is not valid, you will not receive notifications.")

# input
if config['raws']['Metagenomics']:
    MG = [os.path.expandvars(x) if os.path.isabs(os.path.expandvars(x)) else os.getcwd() + "/" + os.path.expandvars(x) for x in config['raws']['Metagenomics'].split()]
else:
    MG = config['raws']['Metagenomics']
if config['raws']['Metatranscriptomics']:
    MT = [os.path.expandvars(x) if os.path.isabs(os.path.expandvars(x)) else os.getcwd() + "/" + os.path.expandvars(x) for x in config['raws']['Metatranscriptomics'].split()]    
else:
    MT = config['raws']['Metatranscriptomics']
if config['raws']['LongReads']:
    if os.path.isabs(os.path.expandvars(config['raws']['LongReads'])):
        LONG = os.path.expandvars(config['raws']['LongReads'])
    else:
        LONG = os.getcwd() + "/" + os.path.expandvars(config['raws']['LongReads'])
else:
    LONG = config['raws']['LongReads']
if config['raws']['Contigs']:
    if os.path.isabs(os.path.expandvars(config['raws']['Contigs'])):
        CONTIGS = os.path.expandvars(config['raws']['Contigs'])
    else:
        CONTIGS = os.getcwd() + "/" + os.path.expandvars(config['raws']['Contigs'])
else:
    CONTIGS = config['raws']['Contigs']
if config['raws']['Alignment_metagenomics']:
    if os.path.isabs(os.path.expandvars(config['raws']['Alignment_metagenomics'])):
        MGaln = os.path.expandvars(config['raws']['Alignment_metagenomics'])
    else:
        MGaln = os.getcwd() + "/" + os.path.expandvars(config['raws']['Alignment_metagenomics'])
else:
    MGaln = config['raws']['Alignment_metagenomics']
if config['raws']['Alignment_metatranscriptomics']:
    if os.path.isabs(os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])):
        MTaln = os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])
    else:
        MTaln = os.getcwd() + "/" + os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])
else:
    MTaln = config['raws']['Alignment_metatranscriptomics']
if config['raws']['Gff']:
    if os.path.isabs(os.path.expandvars(config['raws']['Gff'])):
        ANN = os.path.expandvars(config['raws']['Gff'])
    else:
        ANN = os.getcwd() + "/" + os.path.expandvars(config['raws']['Gff'])
else:
    ANN = config['raws']['Gff']

SAMPLE = config['sample']
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.path.expandvars(config['db_path'])
print(ROOTDIR)
if not os.path.isabs(DBPATH):
    DBPATH = os.path.join(ROOTDIR, DBPATH)
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

# settings
IMP_STEPS = sorted(config['steps'].split())
PREPROCESSING_FILTERING = config['preprocessing_filtering']
FILTER = config['filtering']['filter'].split()
IMP_ASSEMBLER = config['assembly']['assembler']
BIN_STEPS = config['binning']['binners'].split()
if not MG and not MGaln:
    BIN_STEPS = list(set(BIN_STEPS).difference(["MaxBin", "binny"]))
SUM_STEPS = config['summary_steps'].split()

# hardware parameters
BIGCORENO = config['mem']['big_mem_cores']
BIGMEMTOTAL = config['mem']['big_mem_total_gb']
MEMCORE = str(config['mem']['normal_mem_per_core_gb']) + "G"
BIGMEMCORE = str(config['mem']['big_mem_per_core_gb']) + "G"

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = config['tmp_dir']
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

# annotation
# ignore non-required HMMs if annotation was/should be done w/ Mantis
if config["annotation"] == "mantis":
    if "binny" in BIN_STEPS:
        assert "essential" in config["hmm_DBs"].split(), "When using Mantis and Binny you need the \"essential genes\" HMM for annotation"
        config["hmm_DBs"] = "essential"
    else:
        config["hmm_DBs"] = ""
