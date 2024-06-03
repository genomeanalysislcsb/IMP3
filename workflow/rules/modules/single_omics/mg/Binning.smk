# include rules for the Binning part of the workflow

# DAStool
include:
    '../../../Binning/DAStool/dastool.smk'

# MetaWrap
include:
    '../../../Binning/MetaWrap/metawrap.smk'

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']
inputs_binning = ["status/binning_GRiD.done","Binning/MetaWrap/bins.stats"]

## abawaca
if "abawaca" in BIN_STEPS:
    inputs_binning.append("Binning/abawaca/contig2bin.tsv")
    # inputs_binning.append("Binning/abawaca.tar.gz")
    include:
        '../../../Binning/abawaca/abawaca.smk'

## CONCOCT
if "CONCOCT" in BIN_STEPS:
    inputs_binning.append("Binning/CONCOCT/contig2bin.tsv")
    # inputs_binning.append("Binning/CONCOCT.tar.gz")
    include:
        '../../../Binning/CONCOCT/concoct.smk'

## binny
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/contig2bin.tsv")
    # inputs_binning.append("Binning/binny.tar.gz")
    include:
        '../../../Binning/binny/binny_subworkflow.smk'

## MaxBin
if "MaxBin" in BIN_STEPS:
    inputs_binning.append("Binning/MaxBin/maxbin_res.summary")
    # inputs_binning.append("Binning/MaxBin.tar.gz")
    include:
        '../../../Binning/MaxBin/maxbin.smk'
    include:
        '../../../Binning/MaxBin/get_contig2bin.smk'

## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/contig2bin.tsv")
    # inputs_binning.append("Binning/MetaBAT.tar.gz")
    include:
        '../../../Binning/MetaBAT/metabat.smk'

## MetaDecoder
if "MetaDecoder" in BIN_STEPS:
    inputs_binning.append("Binning/MetaDecoder/contig2bin.tsv")
    # inputs_binning.append("Binning/MetaDecoder.tar.gz")
    include:
        '../../../Binning/MetaDecoder/metadecoder.smk'
## SemiBin
if "SemiBin" in BIN_STEPS:
    inputs_binning.append("Binning/SemiBin/contig2bin.tsv")
    # inputs_binning.append("Binning/SemiBin.tar.gz")
    include:
        '../../../Binning/SemiBin/semibin.smk'
## COMEBin
if "COMEBin" in BIN_STEPS:
    inputs_binning.append("Binning/COMEBin/contig2bin.tsv")
    # inputs_binning.append("Binning/COMEBin.tar.gz")
    include:
        '../../../Binning/COMEBin/comebin.smk'

## GRiD
include:
    '../../../Binning/GRiD/grid_dynamic.smk'

if not 'taxonomy' in IMP_STEPS:
    inputs_binning.append("Binning/per_bin_results.tar.gz")


## tar annotation
# if 'analysis' in IMP_STEPS:
#     inputs_binning.append("Analysis/annotation/intermediary.tar.gz")

## master command
rule BINNING:
    input:
        inputs_binning
    output:
        touch('status/binning.done')
