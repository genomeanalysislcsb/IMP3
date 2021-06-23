import os
import shutil
from pathlib import Path


# Snakefile
def open_output(filename):
    return(open(OUTPUTDIR+'/'+filename, 'w+'))

def workflow_ctrl(wildcards):
    return(["status/" + s + ".done" for s in IMP_STEPS if s != "summary"])

def getThreads(max):
    if workflow.cores:
        realThreads = max if max <= workflow.cores else workflow.cores
    elif workflow.nodes:
        realThreads = max if max <= workflow.nodes else workflow.nodes
    else:
        realThreads = max
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
		
# def prepare_input_files(inputs, outputs):
#     """
#     Prepare file names from input into snakemake pipeline.
#     """
#     if len(inputs) != len(outputs):
#         raise OSError("//Inputs and outputs are not of the same length: %s <> %s" % (', '.join(inputs), ', '.join(outputs)))
#     for infilename, outfilename in zip(inputs, outputs):
#         _, fname1 = os.path.split(infilename)
#         _process_file(fname1, infilename, outfilename)
#     # check if r1,r2 are in outputs, if yes, check if both have the same number of lines
#     # quickcheck samtools for bam input
#     # for bam and assembly check number of SQ headers and contigs are equal
#     # for assembly and gff file check that they match

def symlink_input_files(data):
    for infilename, outfilename in data:
        print('symlink',  infilename, '=>', outfilename)
        os.symlink(infilename, outfilename)

def touch_fake_output(files):
    for outfile in files:
        print('touching', outfile)
        Path(outfile).touch()

# def _process_file(fname, inp, outfilename):
#     """
#     Write the input to the output. Handle raw, zip, or bzip input files.
#     """
#     print(inp, '=>', outfilename)
#     import bz2
#     # ungunzip
#     if os.path.splitext(fname)[-1] in ['.gz', '.gzip']:
#         with open(outfilename, 'wb') as whandle, gzip.open(inp, 'rb') as rhandle:
#             whandle.write(rhandle.read())
#     # unbzip2
#     elif os.path.splitext(fname)[-1] in ['.bz2', '.bzip2']:
#         shell("bzip2 -dc {i} > {o}".format(i=inp, o=outfilename))
#     # copy
#     else:
#         shutil.copy(inp, outfilename)

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


def filtering_mg_input(wildcards):
    n = len(wildcards.filterstep.split("."))
    if n == 1:
        return ['Preprocessing/mg.r1.trimmed.fq',
        'Preprocessing/mg.r2.trimmed.fq',
        'Preprocessing/mg.se.trimmed.fq']
    elif n > 1:
        return [s + ".".join([s + "_filtered" for s in FILTER][:(n-1)]) + ".fq" for s in ['Preprocessing/mg.r1.trimmed.','Preprocessing/mg.r2.trimmed.','Preprocessing/mg.se.trimmed.' ]]
    else:
        raise ValueError("invalid filter length %s" % n)

def filtering_mg_SE_input(wildcards):
    n = len(wildcards.filterstep.split("."))
    if n == 1:
        return ['Preprocessing/mg.se.trimmed.fq']
    elif n > 1:
        return ['Preprocessing/mg.se.trimmed.' + ".".join([s + "_filtered" for s in FILTER][:(n-1)]) + ".fq"]
    else:
        raise ValueError("invalid filter length %s" % n)

def filtering_mt_input(wildcards):
    n = len(wildcards.filterstep.split("."))
    if n == 1:
        return ['Preprocessing/mt.r1.trimmed.rna_filtered.fq',
        'Preprocessing/mt.r2.trimmed.rna_filtered.fq',
        'Preprocessing/mt.se.trimmed.rna_filtered.fq']
    elif n > 1:
        return [s + ".".join([s + "_filtered" for s in FILTER][:(n-1)]) + ".fq" for s in ['Preprocessing/mt.r1.trimmed.rna_filtered.','Preprocessing/mt.r2.trimmed.rna_filtered.','Preprocessing/mt.se.trimmed.rna_filtered.' ]]
    else:
        raise ValueError("invalid filter length %s" % n)

def filtering_filter(wildcards):
    return ancient(expand(
            "{dir}/filtering/{filter}.{ext}", filter=wildcards.filterstep.split(".")[-1].split("_")[0],
            ext=['fa', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa'], dir=DBPATH))


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
    if MG or MGaln:
        return expand("Binning/selected_DASTool_bins/{i}/grid",
           i=glob_wildcards("Binning/selected_DASTool_bins/{i}/contigs.fa").i)
    else:
        return ""

def gather_GTDB_bins(wildcards):
    checkpoint_output = checkpoints.getBins.get(**wildcards).output[0]
    return expand("Binning/selected_DASTool_bins/{i}/GTDB",
           i=glob_wildcards("Binning/selected_DASTool_bins/{i}/contigs.fa").i)

def check_bintabs(wildcards):
    checkpoint_output = checkpoints.DAS_tool_check.get(**wildcards).output[0]
    binners = BIN_STEPS
    valid_binners = []
    for binner in binners:
        if os.stat(OUTPUTDIR + "/Binning/" + binner + "/scaffold2bin.tsv").st_size > 0:
            valid_binners.append(binner)
    return expand("Binning/{binner}/scaffold2bin.tsv",binner=valid_binners)
		
def check_bintabs_multi(wildcards):
    checkpoint_output = checkpoints.DASTool_check.get(**wildcards).output[0]
    binners = config['dereplication']['cross_mapping_rebinning']['binners'].split()
    valid_binners = []
    for binner in binners:
        if os.stat(OUTPUTDIR + "/Dereplication/" + binner + "/" + wildcards.sample + "/scaffold2bin.tsv").st_size > 0:
            valid_binners.append(binner)
    return expand("Dereplication/{binner}/{sample}/scaffold2bin.tsv",binner=valid_binners,sample=wildcards.sample)

# MULTI:
def getCatIn(wildcards):
    inputs = ["Catalogue/catalogue_anno.tsv"]
    if config['catalogue']['kraken']['do']:
        inputs.append("Catalogue/catalogue_kraken.parsed.tsv")
    if config['catalogue']['mapping']['do']:
        for type in TYPES:
            inputs.append("Catalogue/" + type + ".counts.tsv")
    return inputs
    
def gather_multi_bins(wildcards):
    checkpoint_output = checkpoints.pre_grid.get(**wildcards).output[0]
    return expand("Dereplication/GRiD_multi/{i}/grid/reads.GRiD.txt",
           i=glob_wildcards("Dereplication/GRiD_multi/{i}/contigs.fa").i)
    
def gather_multi_GTDB_bins(wildcards):
    checkpoint_output = checkpoints.dRep_multi.get(**wildcards).output[0]
    return expand("Dereplication/GTDB_multi/{i}.GTDB_out",
           i=glob_wildcards("Dereplication/dRep_multi_out/dereplicated_genomes/{i}.contigs.fa").i)
           
def getBam(wildcards):
    ori_sample = wildcards.clusterID.split(".")[0]
    return expand("{raw_directory}/Assembly/mg.reads.sorted.bam",
                raw_directory=samples.path[ori_sample])

