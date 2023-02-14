import os
import shutil
from pathlib import Path
import glob
import math

# Snakefile
def open_output(filename):
    return(open(OUTPUTDIR+'/'+filename, 'w+'))

def workflow_ctrl(wildcards):
    return(["status/" + s + ".done" for s in IMP_STEPS if s != "summary"])

def getThreads(cores, bigmem=False):
    if bigmem:
        mem_core = int(config["mem"]["big_mem_per_core"])
        mem_total = int(config["mem"]["big_mem_total"])
        n_cores_from_mem = max(math.floor(mem_total / mem_core), 1)
        n_cores = int(config["mem"]["big_mem_cores"])
        assert n_cores_from_mem == n_cores, f"Number of big_mem cores specified differs from number calculated from " \
                                            f"total memory available and memory per big_mem core: {n_cores} vs " \
                                            f"{n_cores_from_mem}. Please check config.\n" \
                                            f"mem_core: {mem_core}\n" \
                                            f"n_cores: {n_cores}\n" \
                                            f"mem_total: {mem_total}"
        max_cores = int(config["mem"]["big_mem_cores_max"])
    else:
        mem_core = int(config["mem"]["mem_per_core"])
        mem_total = int(config["mem"]["mem_total"])
        n_cores = max(math.floor(mem_total / mem_core), 1)
        max_cores = int(config["mem"]["cores_max"])

    # print(f"mem_core: {mem_core}")
    # print(f"mem_total: {mem_total}")
    # if bigmem:
    #     print(f"n_cores_from_mem: {n_cores_from_mem}")
    # print(f"n_cores: {n_cores}")
    # print(f"max_cores: {max_cores}")
    # print(f"bigmem: {bigmem}")

    # Ensure the maximum number of cores to use per job is actually possible
    assert max_cores <= n_cores, "Maximum number of cores allowed per rule is larger than node core count."

    if workflow.cores:
        realThreads = min(cores, max_cores)
    # What is this for?
    elif workflow.nodes:
        realThreads = min(cores, workflow.nodes)
    else:
        realThreads = max_cores

    # print(f"realThreads: {realThreads}")
    # print("##############################")
    return realThreads

# data input
# inputs (original files) to use
def get_ori_files(wildcards):
    inputs = []
    if MG:
        inputs.append(MG)
    if MT:
        inputs.append(MT)
    if CONTIGS and not "assembly" in IMP_STEPS:
        inputs.append(CONTIGS)
        if MGaln:
            inputs.append(MGaln)
        if MTaln:
            inputs.append(MTaln)
    if "assembly" in IMP_STEPS and LONG:
        inputs.append(LONG)
    if ANN and not "analysis" in IMP_STEPS:
        inputs.append(ANN)
    if any(isinstance(i, list) for i in inputs):
        out = []
        for sublist in inputs:
            if isinstance(sublist, list):
                for item in sublist:
                    out.append(item)
            else:
                out.append(sublist)
    else:
        out = inputs
    return(out)
		
def symlink_input_files(data):
    for infilename, outfilename in data:
        print('symlink',  infilename, '=>', outfilename)
        os.symlink(infilename, outfilename)

def touch_fake_output(files):
    for outfile in files:
        print('touching', outfile)
        Path(outfile).touch()

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
    if wildcards.type in TYPES:
        return input_no_filtering[wildcards.type]
    return 'no_filtering_input_no_file_here'
    
def no_filtering_input_SE(wildcards):
    input_no_filtering = {
        'mg': [
            'Preprocessing/mg.se.trimmed.fq'
        ]
    }
    if wildcards.type == "mg":
        return input_no_filtering[wildcards.type]
    return 'no_filtering_input_no_file_here'

def reads_input_files(wildcards):
    raw_fq_files = [
        'Preprocessing/{type}.r1.fq',
        'Preprocessing/{type}.r2.fq'
    ]
    preprocessed_fq_files = [
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
    ]
    trimmed_fq_files = [
        'Preprocessing/{type}.r1.trimmed.fq',
        'Preprocessing/{type}.r2.trimmed.fq',
        'Preprocessing/{type}.se.trimmed.fq'
    ]
    rna_filtered_fq_mt_files = [
        'Preprocessing/mt.r1.trimmed.rna_filtered.fq',
        'Preprocessing/mt.r2.trimmed.rna_filtered.fq',
        'Preprocessing/mt.se.trimmed.rna_filtered.fq',
    ]
    bam_files = [
        'Assembly/{type}.reads.sorted.bam'
    ]
    if wildcards.type not in TYPES:
        return 'reads_input_files-no-file-here'
    elif MG or MT:
        base = raw_fq_files
        if 'preprocessing' in IMP_STEPS:
            base += preprocessed_fq_files
            if PREPROCESSING_FILTERING:
                base += trimmed_fq_files
                if wildcards.type == 'mt':
                    base += rna_filtered_fq_mt_files
        else:
            if len(MG) == 3 and wildcards.type == 'mg':
                base += ["Preprocessing/mg.se.fq"]
            if len(MT) == 3 and wildcards.type == 'mt':
                base += ["Preprocessing/mt.se.fq"]		
    elif MGaln or MTaln:
        base = bam_files
    return expand(base, type=wildcards.type)
    
def reads_input_files_SE(wildcards):
    raw_fq_files = [
        'Preprocessing/{type}.r1.fq'
    ]
    preprocessed_fq_files = [
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
    ]
    trimmed_fq_files = [
        'Preprocessing/{type}.se.trimmed.fq'
    ]
    if wildcards.type != "mg" :
        return 'reads_input_files-no-file-here'
    else:
        base = raw_fq_files
        base += preprocessed_fq_files
        if PREPROCESSING_FILTERING:
            base += trimmed_fq_files
    return expand(base, type=wildcards.type)

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
    import re
    # list of steps
    filtersteps = wildcards.filterstep.split(".")
    # checks
    assert wildcards.type in TYPES, 'invalid reads type: {} for {}'.format(wildcards.type, wildcards.filterstep)
    if wildcards.type == 'mg':
        assert re.match(r"trimmed", wildcards.filterstep), 'no trimmed in {} for {}'.format(wildcards.filterstep, wildcards.type)
    elif wildcards.type == 'mt':
        assert re.match(r"trimmed\.rna_filtered", wildcards.filterstep), 'no trimmed.rna_filtered in {} for {}'.format(wildcards.filterstep, wildcards.type)
    # input files
    ifiles = expand(
        'Preprocessing/{type}.{read}.{filterstep}.fq',
        type=wildcards.type,
        read=['r1', 'r2', 'se'] if paired else ['se'],
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

# for archiving checkpoint:
def gather_reads(wildcards):
    checkpoint_output = checkpoints.reads_gziped.get(**wildcards).output[0]
    type=wildcards.type
    pot_reads = expand("Preprocessing/{type}.{i}.fq",
           type=type,
           i=glob_wildcards("Preprocessing/"+ type + ".{i}.fq").i)
    last_reads = [os.readlink(l) for l in expand("Preprocessing/{type}.{read}.preprocessed.fq",
                                                   type=type,read=["r1","r2","se"])]
    valid_reads = [f for f in list(set(pot_reads).difference(["Preprocessing/" + s for s in last_reads])) if not os.path.islink(f)]
    return([s + ".gz" for s in valid_reads])


# BINNING:		
# for checkpoints:
def gather_bins(wildcards):
    checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
    # print('gather_bins checkpoint_output:', checkpoint_output)
    if MG or MGaln:
        # print('gather_bins out:', expand("Binning/MetaWrap/bins/{i}/grid", i=glob_wildcards("Binning/MetaWrap/bins/{i}/contigs.fasta").i))
        return expand("Binning/MetaWrap/bins/{i}/grid", i=glob_wildcards("Binning/MetaWrap/bins/{i}/contigs.fasta").i)  # +\
               # expand("Binning/MetaWrap/bins/{i}/grid", i=glob_wildcards("Binning/MetaWrap/bins/{i}/reads.fastq.gz").i)
    else:
        return ""


def gather_GTDB_bins(wildcards):
    checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
    return expand("Binning/MetaWrap/bins/{i}/{tool}", i=glob_wildcards("Binning/MetaWrap/bins/{i}/contigs.fasta").i, tool=["GTDB","grid"])


def check_bintabs(wildcards):
    checkpoint_output = checkpoints.DAS_tool_check.get(**wildcards).output[0]
    binners = BIN_STEPS
    valid_binners = []
    for binner in binners:
        s2b_path = os.path.join(OUTPUTDIR, "Binning", binner, "scaffold2bin.tsv")
        if os.path.exists(s2b_path) and os.stat(s2b_path).st_size > 0:
            valid_binners.append(binner)
    return expand("Binning/{binner}/scaffold2bin.tsv",binner=valid_binners)


def check_binner_bin_fastas(wildcards):
    checkpoint_output = checkpoints.MetaWrap_refine_check.get(**wildcards).output[0]
    # print('check_binner_bin_fastas checkpoint_output:',checkpoint_output)
    binners = BIN_STEPS
    valid_binners = []
    for binner in binners:
        bin_list = glob.glob(OUTPUTDIR + "/Binning/" + binner + "/bins/*.fasta")
        if bin_list:
            if any(os.stat(bin).st_size > 0 for bin in bin_list):
                valid_binners.append(binner)
    if glob.glob(OUTPUTDIR + "/Binning/DASTool/bins/*.fasta"):
        valid_binners.append("DASTool")

    binner_out_paths = expand("Binning/{binner}/bins", binner=valid_binners)
    # print('check_binner_bin_fastas out', binner_out_paths)

    return binner_out_paths


def check_binner_dirs(wildcards):
    checkpoint_output = checkpoints.MetaWrap_refine_check.get(**wildcards).output[0]

    binner_out_paths, = glob_wildcards(os.path.join("Binning", "{binner}"))

    return os.path.join("Binning", "{binner}")



def get_mw_binner_inputs(binner_fasta_dirs, min_comp, max_cont):
    out_list = []
    prefixes = ['-A ', '-B ', '-C ']
    round_1_length = min(len(binner_fasta_dirs), 3)
    input_round_1 = ' '.join([prefix + binner_fasta_path for prefix, binner_fasta_path in
                              zip(prefixes, binner_fasta_dirs[:round_1_length])])
    out_list.append([input_round_1])
    if len(binner_fasta_dirs) > 3:
        input_round_2 = ' '.join([prefix + binner_fasta_path for prefix, binner_fasta_path in
                                  zip(prefixes, [OUTPUTDIR + f'/Binning/MetaWrap/metawrap_round_01/metawrap_{min_comp}_{max_cont}_bins_round_01']
                                                + binner_fasta_dirs[3:5])])
        out_list.append([input_round_2])
    else:
        out_list.append([])
    if len(binner_fasta_dirs) > 5:
        input_round_3 = ' '.join([prefix + binner_fasta_path for prefix, binner_fasta_path in
                                  zip(prefixes, [OUTPUTDIR + f'/Binning/MetaWrap/metawrap_round_02/metawrap_{min_comp}_{max_cont}_bins_round_02']
                                                + binner_fasta_dirs[5:])])
        out_list.append([input_round_3])
    else:
        out_list.append([])
    return out_list


def getBam(wildcards):
    ori_sample = wildcards.clusterID.split(".")[0]
    return expand("{raw_directory}/Assembly/mg.reads.sorted.bam",
                raw_directory=samples.path[ori_sample])


def aggregate_mantis_fc(wildcards):
    checkpoint_output = checkpoints.makegff_mantis.get(**wildcards).output[0]
    return expand("Analysis/annotation/{type}.mantis.{i}.counts.tsv",
           i=glob_wildcards("Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.{i}.gff").i,type=TYPES)


def aggregate_mantis_fc_per_contig(wildcards):
    checkpoint_output = checkpoints.makegff_mantis.get(**wildcards).output[0]
    return expand("Stats/{type}/annotation/{type}.mantis.{i}.counts.tsv",
           i=glob_wildcards("Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.{i}.gff").i,type=TYPES)
