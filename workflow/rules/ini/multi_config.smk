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
import fnmatch



# default configuration file
configfile:
    srcdir("../../../config/config.multi.yaml")


# default executable for snakmake
shell.executable("bash")


# some parameters
SRCDIR = srcdir("../../scripts")
BINDIR = srcdir("../../bin")
ENVDIR = srcdir("../../envs")

# get parameters from the command line
# output
OUTPUTDIR = config['outputdir']

MULTI_STEPS = sorted(config['steps'].split())

try:
    samples = pd.read_csv(os.path.expandvars(config['sample_table']),sep="\t")
except:
    print("Sample table was not found. Please enter the absolute path and file name in the config file.")
    raise
if samples[['sample']].duplicated().any():
    raise Exception('sample names should be unique.')
samples = samples.set_index(["sample"],drop=False)
samples['sample'] = samples['sample'].astype(str)
for sam in samples['sample']:
    if not re.match(r"[0-9a-zA-Z]",sam):
        raise Exception('please start library names with a letter or number')
samples['sample'] = samples['sample'].str.replace('[;|.-]','_')
samples['sample'] = samples['sample'].str.replace('_+','_')
if not 'path' in samples.columns:
    raise Exception("You haven't supplied the paths to the data.")
if samples['path'].duplicated().any():
    raise Exception('paths should be unique.')
if any([not os.path.isdir(path) for path in samples['path'].tolist()]):
    raise Exception('Not all input directories could be found.')
    
multi_config = {}
for sam in samples['sample']:
    spath = samples.loc[samples['sample']==sam, "path"].iat[0] + "/sample.config.yaml"
    if os.path.exists(spath):
#        print("reading configuration from "+ spath)
        with open(spath, 'r') as rhandle:
            multi_config[sam] = yaml.load(rhandle,Loader=yaml.SafeLoader)
    
samples = samples.assign(mg=[1 if (multi_config[sam]['raws']['Metagenomics'].split() or multi_config[sam]['raws']['Alignment_metagenomics'].split() ) else 0 for sam in samples['sample']])
samples = samples.assign(mt=[1 if (multi_config[sam]['raws']['Metatranscriptomics'].split() or multi_config[sam]['raws']['Alignment_metatranscriptomics'].split()) else 0 for sam in samples['sample']])
if "collate" in MULTI_STEPS:
    COLL_STEPS = sorted(config['collection']['results'].split())
    samples = samples.assign(taxonomy=[1 if "taxonomy" in multi_config[sam]['steps'].split() and (multi_config[sam]['raws']['Metagenomics'].split() or multi_config[sam]['raws']['Metatranscriptomics'].split()) else 0 for sam in samples['sample']])
    samples = samples.assign(EukDetect=[1 if "taxonomy" in multi_config[sam]['steps'].split() and (multi_config[sam]['raws']['Metagenomics'].split() or multi_config[sam]['raws']['Metatranscriptomics'].split()) and 'eukdetect' in multi_config[sam] and multi_config[sam]['eukdetect']['run_eukdetect'] else 0 for sam in samples['sample']])
    samples = samples.assign(stats=[1 if "summary" in multi_config[sam]['steps'].split() and "stats" in multi_config[sam]['summary_steps'].split() else 0 for sam in samples['sample']])
    for db in config['hmm_DBs'].split():        
        samples[db] = 0
        samples[db] = [1 if "analysis" in multi_config[sam]['steps'].split() and db in multi_config[sam]['hmm_DBs'].split() else 0 for sam in samples['sample']]
if "catalogue" in MULTI_STEPS:
    samples = samples.assign(genes=[1 if "analysis" in multi_config[sam]['steps'].split() else 0 for sam in samples['sample']])
    samples["prokka_prefix"] = ""
    for sam in samples['sample']:
        if samples.genes[sam] == 1:
            if os.path.isfile(samples.path[sam] + "/Analysis/annotation/prokka.faa"):
                with open(samples.path[sam] + "/Analysis/annotation/prokka.faa",'r') as f:
                    line = f.readline()
                    samples.prokka_prefix[sam] = line.split("_")[0][1:]
if "dereplicate" in MULTI_STEPS:
    samples = samples.assign(bins=[1 if "binning" in multi_config[sam]['steps'].split() else 0 for sam in samples['sample']])
    samples = samples.assign(ass=[1 if ("assembly" in multi_config[sam]['steps'].split() or multi_config[sam]['raws']['Contigs'].split()) else 0 for sam in samples['sample']])
    samples["assembly_type"] = ""
    for sam in samples['sample']:
        if samples.ass[sam] == 1:
            samples.assembly_type[sam] = [f.replace(".assembly.merged.fa","") for f in os.listdir(samples.path[sam] + "/Assembly/") if fnmatch.fnmatch(f,"m*.assembly.merged.fa")][0]
    if config['dereplication']['cross_mapping_rebinning']['do']:
        if not 'rebinning' in samples.columns:
            samples['rebinning'] = "yes"

#print(samples.EukDetect)

TYPES=config['omes'].split()

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))

# temporary directory will be stored inside the OUTPUTDIR directory
# unless an absolute path is set
TMPDIR = os.path.expandvars(config['tmp_dir'])
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

DBPATH = os.path.expandvars(config['db_path'])

# hardware parameters
BIGCORENO = config['mem']['big_mem_cores']
BIGMEMTOTAL = config['mem']['big_mem_total_gb']
MEMCORE = str(config['mem']['normal_mem_per_core_gb']) + "G"
BIGMEMCORE = str(config['mem']['big_mem_per_core_gb']) + "G"

SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

