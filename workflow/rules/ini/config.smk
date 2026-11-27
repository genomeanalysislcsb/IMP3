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

import time
import uuid


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
# Importing steps from the configuration file and sorting them
IMP_STEPS = sorted(config['steps'].split())
# Importing preprocessing_filtering from the configuration file
PREPROCESSING_FILTERING = config['preprocessing_filtering']
FILTER_STEPS = config['filtering']['filter'].split()
# Importing rrna_filtering from the configuration file and ensure type boolean
KEEP_RRNAs = str(config['keep_rrnas']).lower() in ['true', '1', 'yes', 'y', 't']
# Removed options (their conda environments were never installed): stop with a clear message if a config asks for them
_removed = []
if str(config.get("read_coverage_normalization", {}).get("normalize_read_coverage", False)).lower() in ['true', '1', 'yes', 'y', 't']:
    _removed.append("read_coverage_normalization: normalize_read_coverage (BBNorm)")
if str(config.get("eukdetect", {}).get("run_eukdetect", False)).lower() in ['true', '1', 'yes', 'y', 't']:
    _removed.append("eukdetect: run_eukdetect")
if config.get("annotation", "mantis") != "mantis":
    _removed.append("annotation: %s (only mantis)" % config.get("annotation"))
if "DB" in str(config.get("summary_steps", "")).split():
    _removed.append("summary_steps: DB (MongoDB)")
for _b in ("MaxBin", "abawaca"):
    if _b in str(config.get("binning", {}).get("binners", "")).split():
        _removed.append("binning: binners: %s" % _b)
if _removed:
    raise Exception("These IMP3 options were removed, change them in the config: " + "; ".join(_removed))
# Importing filter from the configuration file and splitting it
FILTER = config['filtering']['filter'].split()
# Read filtering against filtering: filter (filtering: method):
#   deacon   (default) every reference with Deacon (minimizers, k 31 w 15). Host genome / transcriptome steps
#            (filtering: deacon_steps) use the human pangenome index (minimizers shared with microbes removed;
#            Constantinides et al. 2025): the first of these steps removes the host reads, later ones pass the reads
#            on. All other references (PhiX, extra filter libraries) use an index built from their fasta
#            (<db>/filtering/deacon/<reference>.k31w15.idx, built on first use)
#   bwa      every reference with bwa-mem2; a read counts as aligned when >= bwa_min_aligned_fraction of it aligns
#            (>= bwa_min_aligned_length bases, mismatches <= bwa_max_mismatch_fraction of the aligned bases); a pair
#            goes when either mate counts (scripts/host_filter_reads.py)
#   bwa_any  IMP3 <= 3.0.3: a pair goes when bwa-mem2 reports either mate as
#            mapped. With short-insert libraries bwa-mem2's pairing / mate rescue marks many microbial pairs as mapped
#            (10 Expobiome stool metagenomes: 12-31 % of the pairs against T2T-CHM13, < 6 % of them human by Deacon
#            and Hostile; bwa and deacon: 0.01-0.1 %)
HOST_METHOD = str(config['filtering'].get('method', 'deacon')).lower()
if HOST_METHOD not in ('deacon', 'bwa', 'bwa_any'):
    raise Exception("filtering: method must be deacon, bwa or bwa_any, not %s" % HOST_METHOD)
DEACON_STEPS = str(config['filtering'].get('deacon_steps', 'T2TCHM13v2genomic T2TCHM13v2transcripts')).split()
DEACON_INDEX = config['filtering'].get('deacon_index', 'deacon/panhuman-1.k31w15.idx')   # under DBPATH/filtering
DEACON_ABS = int(config['filtering'].get('deacon_abs_threshold', 2))
DEACON_REL = float(config['filtering'].get('deacon_rel_threshold', 0.01))
BWA_MIN_FRAC = float(config['filtering'].get('bwa_min_aligned_fraction', 0.8))
BWA_MIN_LEN = int(config['filtering'].get('bwa_min_aligned_length', 30))
BWA_MAX_MM = float(config['filtering'].get('bwa_max_mismatch_fraction', 0.05))
BWA_RULE = "--min-aligned-fraction %s --min-aligned-length %d --max-mismatch-fraction %s" % (BWA_MIN_FRAC, BWA_MIN_LEN,
                                                                                            BWA_MAX_MM)
DEACON_REFS = [s for s in FILTER if s not in DEACON_STEPS]   # with method deacon: own index from the fasta
# Importing assembler from the configuration file
# kraken.memory_mapping: true -> kraken2 --memory-mapping: the database is read through the page cache instead of
# being loaded into the job's memory, and the Kraken2 rules take a whole regular node. Only sensible when the
# database fits into a regular node's memory (otherwise it is re-read from disk all the time). Default false: the
# database is loaded, threads sized from its size and the big-memory specs
KRAKEN_MMAP = str(config.get('kraken', {}).get('memory_mapping', False)).lower() in ['true', '1', 'yes', 'y', 't']
KRAKEN_MMAP_OPT = '--memory-mapping' if KRAKEN_MMAP else ''
IMP_ASSEMBLER = config['assembly']['assembler']
# metaSPAdes (and with it the long-read hybrid assembly) was removed: its hybrid rules never worked, and MEGAHIT
# covers metaG, metaT and the hybrid assembly
if IMP_ASSEMBLER != 'megahit':
    raise Exception("assembly: assembler: %s is not supported; use megahit (metaSPAdes was removed)" % IMP_ASSEMBLER)
if LONG:
    raise Exception("raws: LongReads: long reads were only used by the metaSPAdes assembly, which was removed; "
                    "leave LongReads empty")
# Importing binners from the configuration file and splitting them
BIN_STEPS = config['binning']['binners'].split()
# Checking if the number of binning methods exceeds the maximum limit of 6
if len(BIN_STEPS) > 6:
    raise Exception('A maximum of six binning methods per run are supported at the moment.')
# If neither MG nor MGaln is set, set BIN_STEPS to ['MetaBAT']
if not MG and not MGaln:
    BIN_STEPS = ['MetaBAT']
# Importing summary_steps from the configuration file and splitting them
SUM_STEPS = config['summary_steps'].split()

# hardware parameters
# Defining dictionaries for byte conversion
byte_dict = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
to_gb_dict = {'K': 1e-6, 'M': 1e-3, 'G': 1e0, 'T': 1e3}
to_mb_dict = {'K': 1e-3, 'M': 1e0, 'G': 1e3, 'T': 1e6}

# Importing cluster specifications from the configuration file
MAIN_BIGCORENO = config['cluster_specs']['main_cluster']['big_mem_cores']
# bwa-mem2 reads in batches of 10 Mbp x threads unless -K is given, and its pairing statistics are per batch: the
# filtering and mapping results then depend on the number of threads (a 64-core run removed 2.9 % more read pairs
# as human than the 128-core runs). A fixed batch makes them independent of the threads; the default is the batch
# of a 128-thread run, so results of the usual 128-core runs stay the same
BWA_BATCH = int(config.get('bwa_batch_bases', 1280000000))
MAIN_BIGMEMTOTAL = config['cluster_specs']['main_cluster']['big_mem_total']
MAIN_MEMCORE = str(config['cluster_specs']['main_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_BIGMEMCORE = str(config['cluster_specs']['main_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
MAIN_CORESMAX = config['cluster_specs']['main_cluster']['cores_max']
MAIN_BIGCORESMAX = config['cluster_specs']['main_cluster']['big_mem_cores_max']

# Converting memory units to MB
MAIN_MEMCORE_MB = int(float(MAIN_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
MAIN_BIGMEMCORE_MB = int(float(MAIN_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])

# Define 4 runtime parameter tempalte values to use for all rules (in min)
# 0. very short runtime (2h)
VERY_SHORT_RUNTIME = 120
# 0.1 Shortish runtime (6h)
SHORTISH_RUNTIME = 360
# 1. short runtime (12h)
SHORT_RUNTIME = 720
# 2. medium runtime (24h)
MEDIUM_RUNTIME = 1440
# 3. long runtime (48h)
LONG_RUNTIME = 2880
# 4. very long runtime (120h)
VERY_LONG_RUNTIME = 7200
# 5. ultra long runtime (336h)
ULTRA_LONG_RUNTIME = 20160

# Checking if high_mem_cluster is set in the configuration file
if config['cluster_specs'].get('high_mem_cluster'):
    USE_HIGH_MEM_CLUSTER = True
    # Importing high_mem_cluster specifications from the configuration file
    HIGH_MEM_BIGCORENO = config['cluster_specs']['high_mem_cluster']['big_mem_cores']
    HIGH_MEM_BIGMEMTOTAL = config['cluster_specs']['high_mem_cluster']['big_mem_total']
    HIGH_MEM_MEMCORE = str(config['cluster_specs']['high_mem_cluster']['mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_BIGMEMCORE = str(config['cluster_specs']['high_mem_cluster']['big_mem_per_core']) + config['cluster_specs']['mem_unit']
    HIGH_MEM_CORESMAX = config['cluster_specs']['high_mem_cluster']['cores_max']
    HIGH_MEM_BIGCORESMAX = config['cluster_specs']['high_mem_cluster']['big_mem_cores_max']

    # Converting high_mem_cluster memory units to MB
    HIGH_MEM_MEMCORE_MB = int(float(HIGH_MEM_MEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
    HIGH_MEM_BIGMEMCORE_MB = int(float(HIGH_MEM_BIGMEMCORE[:-1]) * to_mb_dict[config['cluster_specs']['mem_unit']])
else:
    USE_HIGH_MEM_CLUSTER = False

# Get base temporary directory path (usually sample ID)
base_tmpdir = config['tmp_dir']
basename = os.path.basename(base_tmpdir)

# Create unique identifier
unique_id = str(uuid.uuid4())[:8]
unique_tmpdir = f"{basename}_{unique_id}"

# Create full path below the tmp root: /dev/shm (default, RAM), or `tmp_root` in the config / $IMP3_TMP_ROOT for
# machines with little shared memory (workstations, containers)
TMP_ROOT = os.path.expandvars(str(config.get('tmp_root') or os.environ.get('IMP3_TMP_ROOT') or '/dev/shm'))
TMPDIR = os.path.join(TMP_ROOT, unique_tmpdir)

# Create the temporary directory if it doesn't exist
os.makedirs(TMPDIR, exist_ok=True)

# Update all temporary directory references
config['tmp_dir'] = TMPDIR
os.environ['TMPDIR'] = TMPDIR
workflow._shadow_prefix = TMPDIR

# annotation
# If annotation is set to 'mantis' in the configuration file, set hmm_DBs to an empty string
if config['annotation'] == 'mantis':
    config['hmm_DBs'] = ''