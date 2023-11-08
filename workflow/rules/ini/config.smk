import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'

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

if not workflow.config:
    # default configuration file
    configfile:
        srcdir('../../../config/config.imp.yaml')

# default executable for snakmake
shell.executable('bash')

validate(config,schema='../../../schemas/schema.imp.yaml')

# some parameters
SRCDIR = srcdir('../../scripts')
BINDIR = srcdir('../../bin')
ENVDIR = srcdir('../../envs')

# get parameters from the command line
# output
OUTPUTDIR = os.path.abspath(os.path.expandvars(config['outputdir']))
#print(OUTPUTDIR)
SUMDIR = os.path.abspath(os.path.expandvars(config['summarydir']))
if OUTPUTDIR == SUMDIR:
    if SUMDIR[-1] != '/':
        SUMDIR += '_final'
    else:
        SUMDIR = SUMDIR[:-1] + '_final'

EMAIL = config['email']
if EMAIL != '':
    if not re.fullmatch(r'^\w+([\.-]?\w+)*@\w+([\.-]?\w+)*(\.\w{2,3})+$', EMAIL):
        EMAIL = ''
        print('Your email address is not valid, you will not receive notifications.')

# input
if config['raws']['Metagenomics']:
    MG = [os.path.expandvars(x) if os.path.isabs(os.path.expandvars(x)) else os.getcwd() + '/' + os.path.expandvars(x) for x in config['raws']['Metagenomics'].split()]
else:
    MG = config['raws']['Metagenomics']
if config['raws']['Metatranscriptomics']:
    MT = [os.path.expandvars(x) if os.path.isabs(os.path.expandvars(x)) else os.getcwd() + '/' + os.path.expandvars(x) for x in config['raws']['Metatranscriptomics'].split()]
else:
    MT = config['raws']['Metatranscriptomics']
if config['raws']['LongReads']:
    if os.path.isabs(os.path.expandvars(config['raws']['LongReads'])):
        LONG = os.path.expandvars(config['raws']['LongReads'])
    else:
        LONG = os.getcwd() + '/' + os.path.expandvars(config['raws']['LongReads'])
else:
    LONG = config['raws']['LongReads']
if config['raws']['Contigs']:
    if os.path.isabs(os.path.expandvars(config['raws']['Contigs'])):
        CONTIGS = os.path.expandvars(config['raws']['Contigs'])
    else:
        CONTIGS = os.getcwd() + '/' + os.path.expandvars(config['raws']['Contigs'])
else:
    CONTIGS = config['raws']['Contigs']
if config['raws']['Alignment_metagenomics']:
    if os.path.isabs(os.path.expandvars(config['raws']['Alignment_metagenomics'])):
        MGaln = os.path.expandvars(config['raws']['Alignment_metagenomics'])
    else:
        MGaln = os.getcwd() + '/' + os.path.expandvars(config['raws']['Alignment_metagenomics'])
else:
    MGaln = config['raws']['Alignment_metagenomics']
if config['raws']['Alignment_metatranscriptomics']:
    if os.path.isabs(os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])):
        MTaln = os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])
    else:
        MTaln = os.getcwd() + '/' + os.path.expandvars(config['raws']['Alignment_metatranscriptomics'])
else:
    MTaln = config['raws']['Alignment_metatranscriptomics']
if config['raws']['Gff']:
    if os.path.isabs(os.path.expandvars(config['raws']['Gff'])):
        ANN = os.path.expandvars(config['raws']['Gff'])
    else:
        ANN = os.getcwd() + '/' + os.path.expandvars(config['raws']['Gff'])
else:
    ANN = config['raws']['Gff']

SAMPLE = config['sample']
if SAMPLE == '':
    SAMPLE = '_'.join(OUTPUTDIR.split('/')[-2:])
SAMPLE = re.sub('_+','_',re.sub('[;|.-]','_',SAMPLE))
DBPATH = os.path.expandvars(config['db_path'])

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

# settings
IMP_STEPS = sorted(config['steps'].split())
PREPROCESSING_FILTERING = config['preprocessing_filtering']
FILTER = config['filtering']['filter'].split()
IMP_ASSEMBLER = config['assembly']['assembler']
BIN_STEPS = config['binning']['binners'].split()
if len(BIN_STEPS) > 6:
    raise Exception('A maximum of six binning methods per run are supported at the moment.')
if not MG and not MGaln:
    # BIN_STEPS = list(set(BIN_STEPS).difference(['MaxBin', 'binny']))
    BIN_STEPS = ['MetaBAT']
SUM_STEPS = config['summary_steps'].split()

# hardware parameters
byte_dict = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
to_gb_dict = {'K': 1e-6, 'M': 1e-3, 'G': 1e0, 'T': 1e3}
to_mb_dict = {'K': 1e-3, 'M': 1e0, 'G': 1e3, 'T': 1e6}

MAIN_BIGCORENO = config['cluster_specs']['main_cluster']['big_mem_cores']
MAIN_BIGMEMTOTAL = config['cluster_specs']['main_cluster']['big_mem_total']
MAIN_MEMCORE = str(config['cluster_specs']['main_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_BIGMEMCORE = str(config['cluster_specs']['main_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_CORESMAX = config['cluster_specs']['main_cluster']['cores_max']
MAIN_BIGCORESMAX = config['cluster_specs']['main_cluster']['big_mem_cores_max']

MAIN_MEMCORE_MB = int(float(MAIN_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
MAIN_BIGMEMCORE_MB = int(float(MAIN_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])

if config['cluster_specs'].get('high_mem_cluster'):
    print('#' * 10, 'HIGH MEM ON', '#' * 10)
    USE_HIGH_MEM_CLUSTER = True
    HIGH_MEM_BIGCORENO = config['cluster_specs']['high_mem_cluster']['big_mem_cores']
    HIGH_MEM_BIGMEMTOTAL = config['cluster_specs']['high_mem_cluster']['big_mem_total']
    HIGH_MEM_MEMCORE = str(config['cluster_specs']['high_mem_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_BIGMEMCORE = str(config['cluster_specs']['high_mem_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_CORESMAX = config['cluster_specs']['high_mem_cluster']['cores_max']
    HIGH_MEM_BIGCORESMAX = config['cluster_specs']['high_mem_cluster']['big_mem_cores_max']

    HIGH_MEM_MEMCORE_MB = int(float(HIGH_MEM_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
    HIGH_MEM_BIGMEMCORE_MB = int(float(HIGH_MEM_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
else:
    USE_HIGH_MEM_CLUSTER = False

# temporary directory will be stored inside the OUTPUTDIR directory
# unless a absolute path is set
TMPDIR = config['tmp_dir']
if not os.path.isabs(TMPDIR):
    TMPDIR = os.path.join(OUTPUTDIR, TMPDIR)
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

# annotation
# ignore non-required HMMs if annotation was/should be done w/ Mantis
if config['annotation'] == 'mantis':
    config['hmm_DBs'] = ''
