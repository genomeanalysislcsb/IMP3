# Common utils
# See common.rules.smk for common rules

##################################################
# Modules
##################################################
import os
import re
import sys
import math
import yaml
from os import path, getcwd
from os.path import isfile
from itertools import chain

import shutil
import gzip
#import json
import bz2
from copy import deepcopy
import subprocess
import pandas as pd
from pathlib import Path
from snakemake.utils import validate

##################################################
# Variables
##################################################
# TODO

##################################################
# Misc
##################################################

def flatten_list(l):
    """Flatten list (or another iterable) of iterables and return new list with sorted and unique elements"""
    return sorted(set(list(chain.from_iterable(l))))

def add_config_representer():
    """Add YAML representers for config"""
    yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
    yaml.add_representer(tuple, lambda dumper, data: dumper.represent_sequence('tag:yaml.org,2002:seq', data))

def print_config(config_dict=None):
    """Print config in YAML format"""
    if config_dict is None:
        config_dict = config
    print(yaml.dump(config_dict, allow_unicode=True, default_flow_style=False))

def save_config(fname, config_dict=None):
    """Save config to given file"""
    if config_dict is None:
        config_dict = config
    with open(fname, "w") as ofile:
        yaml.dump(config_dict, ofile, allow_unicode=True, default_flow_style=False)

##################################################
# Config
##################################################

def check_input_file_exists(fname):
    """TODO"""
    assert not fname or isfile(fname), "Not an existing regular file: {}".format(fname)

def check_input_file(fname):
    """Input file check"""
    check_input_file_exists(fname)

def check_input_files(config_dict=None):
    """
    check_input_files(config_dict=None)

    Check all input files

    Parameters
    ----------
    config_dict : dict
        Snakemake's config object or a part of it. Set to None to use the complete config.

    Raises
    ------
    Assertion errors if any of the input files does not pass the check
    """
    if config_dict is None:
        config_dict = config["input"].copy()
    for k, v in config_dict.copy().items():
        if isinstance(v, dict):
            check_input_files(v) # recursion
        elif isinstance(v, list):
            for p in v:
                check_input_file(p) # check single file
        else:
            check_input_file(v) # check single file

def get_workflow_default_resources(key):
    """
    workflow_get_default_resources(key=None)

    Get the value of the default resource defined for current the workflow.

    Parameters
    ----------
    key : str

    Returns
    -------
    object or None
        Value saved in default resources or None if the given key does not exist
    """
    return workflow.default_resources.parsed.get(key, None)

def scale_mem(mem, unit="Gb"):
    """Scale memory (in Gb) to Gb, Mb, Kb or B"""
    units = ("Gb", "Mb", "Kb", "B")
    scale = math.pow(1024, units.index(unit)) # will throw an exception if unit not in tuple
    # print("scaling: mem={}, unit={}, scale={} ==> {}".format(mem, unit, scale, mem*scale))
    return int(mem * scale)

def get_mem_big(wildcards, input=None, threads=1, attempt=None, unit="Gb"):
    """
    get_mem_big(wildcards, input=None, threads=1, attempt=None, unit="Gb")

    Callable for snakemake's memory resources for "big memory" nodes
    See: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-resources

    Parameters
    ----------
    wildcards, input, threads, attempt: parameters from snakemake

    unit : str
        Memory unit (Gb, Mg, Kb or B)

    Returns
    -------
    int
        Memory in given units
    """
    return scale_mem(mem=int(get_workflow_default_resources("membig")) * threads, unit=unit)

def get_mem_norm(wildcards, input=None, threads=1, attempt=None, unit="Gb"):
    """
    get_mem_norm(wildcards, input=None, threads=1, attempt=None, unit="Gb")

    Callable for snakemake's memory resources for "normal memory" nodes
    See: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-resources

    Parameters
    ----------
    wildcards, input, threads, attempt: parameters from snakemake

    unit : str
        Memory unit (Gb, Mg, Kb or B)

    Returns
    -------
    int
        Memory in given units
    """
    return scale_mem(mem=int(get_workflow_default_resources("memnorm")) * threads, unit=unit)

def get_cores_big():
    """Get number of big-mem cores"""
    return int(get_workflow_default_resources("coresbig"))

def get_cores_norm():
    """Get number of big-mem cores"""
    return int(get_workflow_default_resources("coresnorm"))

def get_abs_path(pathname):
    """
    get_abs_path(pathname)

    Get the absolute path from given (relative) path.
    Given path can also include environment variables which will be expanded.
    If the path is an empty string or None, None will be returned.

    Parameters
    ----------
    pathname : str

    Returns
    -------
    str or None
        Absolute path or None if no/empty path was given.
    """
    # type check
    assert pathname is None or isinstance(pathname, str), "Path name is not None or string: {}".format(pathname)
    if not pathname: # TODO: return None
        return ""
        # return None
    # expand
    pathname = path.expandvars(pathname)
    if not path.isabs(pathname):
        # print("get_abs_path({}): CWD = {}".format(pathname, getcwd()))
        # pathname = path.join(getcwd(), pathname)
        pathname = path.abspath(pathname)
    return pathname

def create_sample_name(pathname):
    """
    create_sample_name(pathname)

    Derive sample name from given path (supposed to be the output/working directory).
    The last two directories will be joined by an underscore.
    Any non-word characters, i.e. those not in [a-zA-Z0-9_], will be removed.

    Parameters
    ----------
    pathname : str

    Returns
    -------
    str
        Sample name
    """
    # derive from given path
    sample = "_".join(os.path.normpath(path).split(os.sep)[-2:])
    # string cleaning
    sample = re.sub(r"\W", "", sample)
    return sample

def get_steps():
    """Get workflow steps"""
    return config["steps"]

def proc_input_file(fname):
    """
    TODO
    """
    fname = get_abs_path(fname)
    return fname

def proc_input_files(config_dict=None):
    """
    TODO
    """
    if config_dict is None:
        config_dict = config["input"].copy()
    # print("\n\nProcess\n{}".format(print_config(config_dict)))
    for k, v in config_dict.copy().items():
        if isinstance(v, dict):
            config_dict[k] = proc_input_files(v)
        elif isinstance(v, list):
            config_dict[k] = [proc_input_file(p) for p in v]
        else:
            config_dict[k] = proc_input_file(v)
    return config_dict

def has_sr(omic, reads=None):
    """XXX"""
    if reads is None:
        return any(config["input"]["sr"][omic].values())
    elif reads == "pe":
        return bool(config["input"]["sr"][omic]["r1"]) and bool(config["input"]["sr"][omic]["r2"])
    elif reads == "se":
        return bool(config["input"]["sr"][omic]["se"])

def has_lr(omic):
    """XXX"""
    return omic in config["input"]["lr"].keys() and bool(config["input"]["lr"][omic])

def has_asm_map(omic):
    """XXX"""
    return bool(config["input"]["assembly"]["mapping"][omic])

def has_asm():
    """XXX"""
    return bool(config["input"]["assembly"]["contigs"])

def has_asm_gff():
    """XXX"""
    return bool(config["input"]["assembly"]["gff"])

def get_omics():
    """
    get_omics()

    Get omics types based on user input data.

    Returns
    -------
    tuple of str
        Tuple of omics types (strings; possible values = {mg, mt})
    """
    if has_sr("mg") or has_sr("mt"):
        return tuple([k for k in ["mg", "mt"] if has_sr(k)])
    if has_asm_map("mg") or has_asm_map("mt"):
        return tuple([k for k in ["mg", "mt"] if has_asm_map(k)])
    assert True, "Cannot determine the omics type: need either reads and/or assembly mapping as input"
    # if has_asm() and (has_sr("mg") or has_sr("mt")):
    #     return tuple([k for k in ["mg", "mt"] if has_sr(k)])

def has_mg():
    """Has MG data"""
    return "mg" in get_omics()

def has_mt():
    """Has MT data"""
    return "mt" in get_omics()

def has_mgmt():
    """Has MG/MT data"""
    return "mg" in get_omics() and "mt" in get_omics()

def get_omics_assembly():
    """
    get_omics_assembly()

    Get assembly type based on user input data.

    Returns
    -------
    str
        Assembly type (possible values = {mgmt, mg, mt})
    """
    # NOTE: not checking for PE/SE !!!
    if has_sr("mg") and has_sr("mt") and config["assembly"]["hybrid"]:
        return "mgmt"
    if has_sr("mg"): # even if MT is also given
        return "mg"
    if has_sr("mt"):
        return "mt"
    if has_asm_map("mg"): # even if MT is also given
        return "mg"
    if has_asm_map("mt"):
        return "mt"
    if has_asm() and has_sr("mg"): # even if MT is also given
        return "mg"
    if has_asm() and has_sr("mt"):
        return "mt"

def get_read_types_raw_(omic):
    """TODO"""
    if not has_sr(omic, "pe"): # SE only
        return ["se"]
    else:
        return ["r1", "r2"]

def get_read_types_raw(omic=None):
    """TODO"""
    omic = get_omics() if omic is None else [omic]
    return flatten_list([get_read_types_raw_(omic) for omic in get_omics()])

def get_read_types_proc_(omic):
    """TODO"""
    if not has_sr(omic, "pe"): # SE only
        return ["se"]
    else:
        return ["r1", "r2", "se"]

def get_read_types_proc(omic=None):
    """TODO"""
    omic = get_omics() if omic is None else [omic]
    return flatten_list([get_read_types_proc_(omic) for omic in get_omics()])

def get_read_types_trimmoed_(omic):
    """TODO"""
    if not has_sr(omic, "pe"): # SE only
        return ["se"]
    else: # PE + SE
        return ["r1", "r2", "se1", "se2"]

def get_read_types_trimmoed(omic=None):
    """TODO"""
    omic = get_omics() if omic is None else [omic]
    return flatten_list([get_read_types_trimmoed_(omic) for omic in get_omics()])

def get_read_types_trimmed_(omic):
    """TODO"""
    reads = ["r1", "r2", "se1", "se2", "se"]
    if not has_sr(omic, "pe"): # SE only
        reads = ["se"]
    elif config["preprocessing"]["nextseq"]: # PE + SE
        reads += ["r1Cut", "r2Cut", "se1Cut", "se2Cut"]
    return reads

def get_read_types_trimmed(omic=None):
    """TODO"""
    omic = get_omics() if omic is None else [omic]
    return flatten_list([get_read_types_trimmed_(omic) for omic in get_omics()])

def get_filtering_genomes():
    """Get genomes to be used for reads filtering"""
    return config["preprocessing"]["filtering"]["genomes"]

def do_filtering_genomes():
    """Whether any genomes were given to filter reads"""
    return len(get_filtering_genomes()) > 0

def get_assembler():
    return config["assembly"]["assembler"]

# TODO: Q: Why should this function be used? Snakemake is supposed to take care of the number of threads to be used.
# def getThreads(max):
#     if workflow.cores:
#         realThreads = max if max <= workflow.cores else workflow.cores
#     elif workflow.nodes:
#         realThreads = max if max <= workflow.nodes else workflow.nodes
#     else:
#         realThreads = max
#     return realThreads

##################################################
# Functions for snakemake rules
##################################################

# def workflow_ctrl(wildcards):
#     return(["status/" + s + ".done" for s in get_steps() if s != "summary"])

def get_ori_files(config_dict=None):
    """Get all input data files"""
    files = []
    if config_dict is None:
        config_dict = config["input"].copy()
    for k, v in config_dict.items():
        if isinstance(v, dict):
            files += get_ori_files(v)
        elif isinstance(v, list):
            files += [p for p in v if p]
        elif v:
            files.append(v)
    return files

def get_ori_files_targets():
    """Get targets for all input data files"""
    files = []
    # reads
    if has_sr("mg") or has_sr("mt"):
        for omic in config["input"]["sr"].keys():
            for reads, fname in config["input"]["sr"][omic].items():
                if fname:
                    files.append("Preprocessing/{omic}.{reads}.fq".format(omic=omic, reads=reads))
    # long reads
    if "assembly" in get_steps() and has_lr("mg"):
        files.append("Assembly/long_reads.fastq")
    # assembly + mapping
    if "assembly" not in get_steps() and has_asm():
        files.append("Assembly/%.assembly.merged.fa" % get_omics_assembly())
        for omic, fname in config["input"]["assembly"]["mapping"].items():
            if fname:
                files.append("Assembly/{omic}.reads.sorted.bam".format(omic=omic))
    if "analysis" not in get_steps() and has_asm_gff():
        files.append("Analysis/annotation/annotation_CDS_RNA_hmms.gff")
    return files

def get_ori_files_symlinks_sr_se(omic, reads):
    """TODO"""
    if not has_sr(omic, reads):
        return "Preprocessing/{omic}.{reads}.fq".format(omic=omic, reads=reads)

def get_ori_files_symlinks():
    """TODO"""
    files = []
    files += get_ori_files_symlinks_sr_se()
    return files

def get_read_files_raw(omic):
    """Raw FASTQ files (MG or MT)"""
    return expand(
        "Preprocessing/{omic}.{reads}.fq",
        omic=omic, reads=get_read_types_raw(omic),
    )

def get_read_files_trimmed(omic, only_final=False):
    """Trimmed FASTQ files (MG or MT)"""
    if only_final:
        return expand(
            "Preprocessing/{omic}.{reads}.trimmed.fq",
            omic=omic, reads=get_read_types_proc(omic)
        )
    files = expand(
        "Preprocessing/{omic}.{reads}.trimmed.fq",
        omic=omic, reads=get_read_types_trimmed(omic)
    )
    # for nextseq only
    if config["preprocessing"]["nextseq"]:
        files += expand(
            "Preprocessing/{omic}.{reads}.trimmoed.fq",
            omic=omic, reads=get_read_types_trimmoed(omic)
        )
    return files

def get_read_files_filtered_rrna():
    """Read FASTQ files (MT) after rRNA filtering"""
    return expand(
        "Preprocessing/{omic}.{reads}.trimmed.rna_filtered.fq",
        omic="mt", reads=get_read_types_proc("mt")
    )

def get_read_files_filtered_genomes_basename(i=None, use_re=False):
    """Basename for a filtered FASTQ file"""
    genomes = get_filtering_genomes()
    if i is not None:
        genomes = genomes[:i]
    sep="."
    if use_re:
        sep="\."
    return sep.join([g + "_filtered" for g in genomes])

def get_read_files_filtered_genomes_basenames(use_re=False):
    """Basenames for filtered FASTQ files"""
    return [get_read_files_filtered_genomes_basename(x, use_re) for x in range(1, len(get_filtering_genomes()) + 1)]

def get_read_files_filtered_genomes(omic):
    """FASTQ files (MG or MT) after genome-filtering"""
    return expand(
        "Preprocessing/{omic}.{reads}.{procstep}.{filtsteps}.fq",
        omic=omic, reads=get_read_types_proc(omic),
        procstep="trimmed.rna_filtered" if omic == "mt" else "trimmed",
        filtsteps=get_read_files_filtered_genomes_basenames(False)
    )

def get_read_files_proc(omic):
    """Final preprocessed FASTQ files (MG or MT)"""
    return expand(
        "Preprocessing/{omic}.{reads}.preprocessed.fq",
        omic=omic, reads=get_read_types_proc(omic),
    )

def get_read_files_all(exclude_proc=True):
    """Get all (to be) generated FASTQ files"""
    files = []
    for omic in get_omics():
        files += get_read_files_raw(omic)
        files += get_read_files_trimmed(omic)
        if omic == "mt":
            files += get_read_files_filtered_rrna()
        files += get_read_files_filtered_genomes(omic)
        if not exclude_proc:
            files += get_read_files_proc(omic)
    return files

# # def get_read_files(omic=None, all_steps=False):
# #     """
# #     TODO
# #     """
# #     if omic is None:
# #         omic = get_omics()
# #     # raw
# #     files = expand(
# #         "Preprocessing/{omic}.{reads}.fq",
# #         omic=omic, reads=get_read_types(omic)
# #     )
# #     # processed files
# #     if XXX:
# #         files += expand(
# #             "Preprocessing/{omic}.{reads}.processed.fq",
# #             omic=omic, reads=get_read_types(omic)
# #         )

# #     # processed (trimmed, RNA-filtered)
# #     rtypes = []
# #     if "preprocessing" in get_steps():
# #         rtypes.append("preprocessed")
# #         if "mt" in omic:
# #             rtypes.append("trimmed.rna_filtered")
# #         if do_filtering_genomes():
# #             rtypes.append("trimmed")
# #             if all_steps:
# #                 if "mg" in omic:
# #                     rtypes += get_read_files_filtered_genomes_basenames("", False)
# #                 if "mt" in omic:
# #                     rtypes += get_read_files_filtered_genomes_basenames("", False)
# #         files += expand(
# #                 "Preprocessing/{omic}.{reads}.{rtype}.fq",
# #                 omic=omic, reads=get_read_types(omic), rtype=rtypes
# #         )
# #     return files

def get_read_files_wrapper(wildcards):
    """Wrapper for snakemake rules"""
    return get_read_files(wildcards.omic)

def get_trimmomatic_params(wildcards):
    """Get Trimmomatic parameters for a snakemake rule"""
    return "ILLUMINACLIP:{adapter}:{smm}:{pct}:{sct} LEADING:{lead} TRAILING:{trail} MINLEN:{minlen} SLIDINGWINDOW:{ws}:{wq} MAXINFO:{tl}:{strict}".format(
        adapter=os.path.join(DBPATH, "adapters", "%s.fa" % config["preprocessing"]["trimmomatic"]["adapter"][wildcards.omic]),
        smm=config["preprocessing"]["trimmomatic"]["seed_mismatch"],
        pct=config["preprocessing"]["trimmomatic"]["palindrome_clip_threshold"],
        sct=config["preprocessing"]["trimmomatic"]["simple_clip_threshold"],
        lead=config["preprocessing"]["trimmomatic"]["leading"],
        trail=config["preprocessing"]["trimmomatic"]["trailing"],
        minlen=config["preprocessing"]["trimmomatic"]["minlen"],
        ws=config["preprocessing"]["trimmomatic"]["window_size"],
        wq=config["preprocessing"]["trimmomatic"]["window_quality"],
        strict=config["preprocessing"]["trimmomatic"]["strictness"],
        tl=config["preprocessing"]["trimmomatic"]["target_length"],
    )

def get_cutadapt_nextset_params(wildcards):
    """Get parameters for NextSeq trimming w/ Cutadapt"""
    return "--nextseq-trim={} --minimum-length {}".format(
        config["preprocessing"]["trimmomatic"]["trailing"],
        config["preprocessing"]["trimmomatic"]["minlen"]
    )

# PREPROCESSING:
def no_filtering_input(wildcards):
    input_no_filtering = {
        'mg': [
            'Preprocessing/mg.r1.trimmed.fq',
            'Preprocessing/mg.r2.trimmed.fq',
            'Preprocessing/mg.se.trimmed.fq'
        ],
        'mt': [
            'Preprocessing/mt.r1.trimmed.rna_filtered.fq',
            'Preprocessing/mt.r2.trimmed.rna_filtered.fq',
            'Preprocessing/mt.se.trimmed.rna_filtered.fq'
        ]
    }
    if wildcards.omic in get_omics():
        return input_no_filtering[wildcards.omic]
    return 'no_filtering_input_no_file_here'
    
def no_filtering_input_SE(wildcards):
    input_no_filtering = {
        'mg': [
            'Preprocessing/mg.se.trimmed.fq'
        ]
    }
    if wildcards.omic == "mg":
        return input_no_filtering[wildcards.omic]
    return 'no_filtering_input_no_file_here'

def filtering_reads_pe_input(wildcards):
    """Input files for paired-end MG/MT reads filtering (excl. rRNA filtering)"""
    return filtering_reads_input(wildcards, True)

def filtering_reads_se_input(wildcards):
    """Input files for single-end MG/MT reads filtering (excl. rRNA filtering)"""
    return filtering_reads_input(wildcards, False)

def filtering_reads_input(wildcards, paired=True):
    """
    FASTQ input files for the paired-/single-end MG/MT reads filtering rule.
    Note, that rRNA filtering for MT reads has a separate rule.

    Parameters
    ----------
    wildcards :
        wildcards from the rule, must contain keywords 'type', 'filterstep'
    paired : bool
        paired-end (R1, R2, SE) or single-end (SE) reads
    """
    # list of steps
    filtersteps = wildcards.filterstep.split(".")
    # checks
    assert wildcards.omic in get_omics(), 'invalid reads type: {} for {}'.format(wildcards.omic, wildcards.filterstep)
    if wildcards.omic == 'mg':
        assert re.match(r"trimmed", wildcards.filterstep), 'no trimmed in {} for {}'.format(wildcards.filterstep, wildcards.omic)
    elif wildcards.omic == 'mt':
        assert re.match(r"trimmed\.rna_filtered", wildcards.filterstep), 'no trimmed.rna_filtered in {} for {}'.format(wildcards.filterstep, wildcards.omic)
    # input files
    ifiles = expand(
        'Preprocessing/{omic}.{reads}.{filterstep}.fq',
        omic=wildcards.omic,
        reads=['r1', 'r2', 'se'] if paired else ['se'],
        filterstep='.'.join(filtersteps[:-1])
    )
    # print("{}: {}".format(filtersteps, ifiles))
    return ifiles

def filtering_filter(wildcards):
    """Reference input files for the paired-/single-end MG/MT reads filtering rule."""
    return ancient(expand(
        "{dir}/filtering/{filter}.{ext}",
        dir=DBPATH,
        filter=wildcards.filterstep.split(".")[-1].split("_")[0],
        ext=['fa', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa']
    ))

# # for archiving checkpoint:
# def gather_reads(wildcards):
#     checkpoint_output = checkpoints.reads_gziped.get(**wildcards).output[0]
#     omic=wildcards.omic
#     pot_reads = expand("Preprocessing/{omic}.{i}.fq",
#            omic=omic,
#            i=glob_wildcards("Preprocessing/"+ omic + ".{i}.fq").i)
#     last_reads = [os.readlink(l) for l in expand("Preprocessing/{omic}.{reads}.preprocessed.fq",
#                                                    omic=omic,read=["r1","r2","se"])]
#     valid_reads = [f for f in list(set(pot_reads).difference(["Preprocessing/" + s for s in last_reads])) if not os.path.islink(f)]
#     return([s + ".gz" for s in valid_reads])

# # BINNING:		
# # for checkpoints:
# def gather_bins(wildcards):
#     checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
#     if MG or MGaln:
#         return expand("Binning/selected_DASTool_bins/{i}/grid",
#            i=glob_wildcards("Binning/selected_DASTool_bins/{i}/contigs.fa").i)
#     else:
#         return ""

# def gather_GTDB_bins(wildcards):
#     checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
#     return expand("Binning/selected_DASTool_bins/{i}/GTDB",
#            i=glob_wildcards("Binning/selected_DASTool_bins/{i}/contigs.fa").i)

# def check_bintabs(wildcards):
#     checkpoint_output = checkpoints.DAS_tool_check.get(**wildcards).output[0]
#     binners = BIN_STEPS
#     valid_binners = []
#     for binner in binners:
#         if os.stat(OUTPUTDIR + "/Binning/" + binner + "/scaffold2bin.tsv").st_size > 0:
#             valid_binners.append(binner)
#     return expand("Binning/{binner}/scaffold2bin.tsv",binner=valid_binners)
		
# def check_bintabs_multi(wildcards):
#     checkpoint_output = checkpoints.DASTool_check.get(**wildcards).output[0]
#     binners = config['dereplication']['cross_mapping_rebinning']['binners'].split()
#     valid_binners = []
#     for binner in binners:
#         if os.stat(OUTPUTDIR + "/Dereplication/" + binner + "/" + wildcards.sample + "/scaffold2bin.tsv").st_size > 0:
#             valid_binners.append(binner)
#     return expand("Dereplication/{binner}/{sample}/scaffold2bin.tsv",binner=valid_binners,sample=wildcards.sample)

# # MULTI:
# def getCatIn(wildcards):
#     inputs = ["Catalogue/catalogue_anno.tsv"]
#     if config['catalogue']['kraken']['do']:
#         inputs.append("Catalogue/catalogue_kraken.parsed.tsv")
#     if config['catalogue']['mapping']['do']:
#         for omic in get_omics():
#             inputs.append("Catalogue/" + omic + ".counts.tsv")
#     return inputs
    
# def gather_multi_bins(wildcards):
#     checkpoint_output = checkpoints.pre_grid.get(**wildcards).output[0]
#     return expand("Dereplication/GRiD_multi/{i}/grid/reads.GRiD.txt",
#            i=glob_wildcards("Dereplication/GRiD_multi/{i}/contigs.fa").i)
    
# def gather_multi_GTDB_bins(wildcards):
#     checkpoint_output = checkpoints.dRep_multi.get(**wildcards).output[0]
#     return expand("Dereplication/GTDB_multi/{i}.GTDB_out",
#            i=glob_wildcards("Dereplication/dRep_multi_out/dereplicated_genomes/{i}.contigs.fa").i)
           
# def getBam(wildcards):
#     ori_sample = wildcards.clusterID.split(".")[0]
#     return expand("{raw_directory}/Assembly/mg.reads.sorted.bam",
#                 raw_directory=samples.path[ori_sample])

def get_bwa_idx_ext():
    """BWA index file extensions"""
    return ["amb", "ann", "bwt", "pac", "sa"]

def get_prokka_ext(which="all"):
    """Prokka output file extensions (all, FASTA, ...)"""
    # https://github.com/tseemann/prokka#output-files
    if which == "fasta":
        return ["faa", "fna", "ffn", "fsa"]
    return ["gff", "gbk", "fna", "faa", "ffn", "sqn", "fsa", "tbl", "err", "log", "txt", "tsv"]

def get_hmms(get_path=False):
    """Get HMM names or paths"""
    hmms = config["analysis"]["hmm"]
    if get_path:
        hmms = [get_hmm_path(hmm) for hmm in hmms]
    return hmms

def get_hmm_path(wildcards):
    """Get path to HMM database file"""
    return os.path.join(DBPATH, "hmm", wildcards.db)

def get_hmm_cutoff(wildcards):
    """Get hmm cutoff parameters"""
    return config["analysis"]["hmm_settings"][wildcards.db]["cutoff"]

def get_hmm_trim(wildcards):
    """Get hmm trimming parameters"""
    return config["analysis"]["hmm_settings"][wildcards.db]["trim"]

def get_feature_counts_stranded(omic):
    """Get the 'stranded' parameter for featureCounts"""
    return config["analysis"]["feature_counts"]["stranded"][omic]

def get_rrna_refs():
    """Get rRNA reference files"""
    return expand("{dbpath}/sortmerna/{files}.fasta", dbpath=DBPATH, files=config["preprocessing"]["filtering"]["rrna"])

def get_rrna_refs_idx():
    """Get rRNA reference files sortmerna index files"""
    return expand(
        "{dbpath}/sortmerna/idx/{files}.{ext}",
        dbpath=DBPATH,
        files=config["preprocessing"]["filtering"]["rrna"],
        ext=["bursttrie_0.dat", "kmer_0.dat", "pos_0.dat", "stats"]
    )

def get_mem_samtools_sort(threads=1):
    """Get mem for samtools sort"""
    # mem * cores * 0.75 - 0.5
    return get_mem_big(threads) * 0.75 - 0.5

def get_mapping_header(wildcards):
    return "@RG\\tID:{sample}\\tSM:{omic}".format(sample=config["sample"], omic=wildcards.omic)

def metaspades_lr_input(wildcards):
    """XXX"""
    if not has_lr(wildcards.omic):
        return []
    else:
        return ["Assembly/long_reads.fastq"]

def metaspades_lr_input_params(wildcards):
    """XXX"""
    if not has_lr(wildcards.omic):
        return ""
    else:
        return "--{lr_tech}".format(lr_tech=config["lr_tech"])

def get_host_proteome():
    """XXX"""
    return config["proteomics"]["host_proteome"]

def has_host_proteome():
    """XXX"""
    return bool(get_host_proteome())