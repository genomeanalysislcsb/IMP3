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
    srcdir("../../../config/config.imp.yaml")


# default executable for snakmake
shell.executable("bash")

# custom configuration file
CUSTOM_CONFIG_PATH = os.environ.get("CONFIGFILE", "../../../config/userconf.imp.yaml")

# some relevant paths
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")
CONFIGDIR = srcdir("../../../config")
MODULDIR = srcdir("../../../submodules")

# merge 2 configurations files together
if os.path.exists(CUSTOM_CONFIG_PATH):
    #print("reading configuration from "+CUSTOM_CONFIG_PATH)
    with open(CUSTOM_CONFIG_PATH, 'r') as rhandle:
        data = yaml.load(rhandle)
        config = dict_merge(config, data)

# PathoFact -- BEGIN
# Relevant for input/output data
config["pathofact"]["sample"] = ["pathofact"]
config["pathofact"]["project"] = "pathofact"
config["pathofact"]["workflow"] = "complete"
config["pathofact"]["datadir"] = os.path.join(
    os.path.abspath(os.path.expandvars(os.environ.get("OUTPUTDIR", config['outputdir']))),
    "Analysis", "annotation"
)
# Paths for data/tools
config["pathofact"]["scripts"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/scripts"))
config["pathofact"]["deeparg"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/submodules/deeparg-ss/deepARG.py"))
config["pathofact"]["deepvirfinder"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/submodules/DeepVirFinder/dvf.py"))
config["pathofact"]["tox_hmm"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/databases/toxins/combined_Toxin.hmm"))
config["pathofact"]["tox_lib"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/databases/library_HMM_Toxins.csv"))
config["pathofact"]["vir_hmm"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/databases/virulence/Virulence_factor.hmm"))
config["pathofact"]["vir_domains"] = os.path.abspath(os.path.join(MODULDIR, "PathoFact/databases/models_and_domains"))
# Set mem
config["pathofact"]["mem"] = dict()
config["pathofact"]["mem"]["normal_mem_per_core_gb"] = str(os.environ.get("MEMCORE", config['mem']['normal_mem_per_core_gb'])) + "G"
config["pathofact"]["mem"]["big_mem_cores"] = os.environ.get("BIGCORENO", config['mem']['big_mem_cores'])
config["pathofact"]["mem"]["big_mem_per_core_gb"] = str(os.environ.get("BIGMEMCORE", config['mem']['big_mem_per_core_gb'])) + "G"
# Parameters for PlasFlow
config["pathofact"]["plasflow_threshold"] = config["plasflow"]["threshold"]
config["pathofact"]["plasflow_minlen"] = config["plasflow"]["minlen"]
# print(config["pathofact"])
# PathoFact -- END

# config validation
validate(config,schema="../../../schemas/schema.imp.yaml")

# some parameters
# SRCDIR = srcdir("../../scripts")
# BINDIR = srcdir("../../bin")
# ENVDIR = srcdir("../../envs")
# CONFIGDIR = srcdir("../../../config")

# get parameters from the command line
# output
OUTPUTDIR = os.path.abspath(os.path.expandvars(os.environ.get("OUTPUTDIR", config['outputdir'])))
#print(OUTPUTDIR)
SUMDIR = os.path.abspath(os.path.expandvars(os.environ.get("SUMDIR", config['summarydir'])))
if OUTPUTDIR == SUMDIR:
    if SUMDIR[-1] != "/":
        SUMDIR += "_final"
    else:
        SUMDIR = SUMDIR[:-1] + "_final"

# input
MG = os.environ.get("MG", config['raws']['Metagenomics']).split()
MT = os.environ.get("MT", config['raws']['Metatranscriptomics']).split()
LONG = os.environ.get("LONGREADS", config['raws']['LongReads'])
CONTIGS = os.environ.get("CONTIGS", config['raws']['Contigs'])
MGaln = os.environ.get("MGaln", config['raws']['Alignment_metagenomics'])
MTaln = os.environ.get("MTaln", config['raws']['Alignment_metatranscriptomics'])
ANN = os.environ.get("ANN", config['raws']['Gff'])
SAMPLE = os.environ.get("SAMPLE", config['sample'])
if SAMPLE == "":
    SAMPLE = "_".join(OUTPUTDIR.split("/")[-2:])
SAMPLE = re.sub("_+","_",re.sub("[;|.-]","_",SAMPLE))
DBPATH = os.environ.get("DBPATH", os.path.expandvars(config['db_path']))
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
IMP_STEPS = sorted(os.environ.get("IMP_STEPS", config['steps']).split())
PREPROCESSING_FILTERING = os.environ.get("PREPROCESSING_FILTERING", config['preprocessing_filtering'])
FILTER = os.environ.get("FILTER", config['filtering']['filter']).split()
IMP_ASSEMBLER = os.environ.get("IMP_ASSEMBLER", config['assembly']['assembler'])
BIN_STEPS = os.environ.get("IMP_BINNING_METHOD", config['binning']['binners']).split()
if not MG and not MGaln:
    BIN_STEPS = list(set(BIN_STEPS).difference(["MaxBin", "binny"]))
SUM_STEPS = os.environ.get("SUM_STEPS", config['summary_steps']).split()

# hardware parameters
BIGCORENO = os.environ.get("BIGCORENO", config['mem']['big_mem_cores'])
BIGMEMTOTAL = os.environ.get("BIGMEMTOTAL", config['mem']['big_mem_total_gb'])
MEMCORE = str(os.environ.get("MEMCORE", config['mem']['normal_mem_per_core_gb'])) + "G"
BIGMEMCORE = str(os.environ.get("BIGMEMCORE", config['mem']['big_mem_per_core_gb'])) + "G"

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = os.environ.get("TMPDIR", config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)
