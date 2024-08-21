# include rules for the Binning part of the workflow

include:
    '../../Binning/DAStool/dastool.smk'

# MetaWrap
include:
    '../../Binning/MetaWrap/metawrap.smk'

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']

inputs_binning = ["status/binning_GRiD.done", "Binning/MetaWrap/bins.stats", "Binning/DASTool"]

## abawaca
if "abawaca" in BIN_STEPS:
    inputs_binning.append("Binning/abawaca/contig2bin.tsv")
    include:
        '../../Binning/abawaca/abawaca.smk'

## CONCOCT
if "CONCOCT" in BIN_STEPS:
    inputs_binning.append("Binning/CONCOCT/contig2bin.tsv")
    include:
        '../../Binning/CONCOCT/concoct.smk'

## binny
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/contig2bin.tsv")
    include:
        '../../Binning/binny/binny_subworkflow.smk'

## MaxBin
if "MaxBin" in BIN_STEPS:
    inputs_binning.append("Binning/MaxBin/maxbin_res.summary")
    include:
        '../../Binning/MaxBin/maxbin.smk'
    include:
        '../../Binning/MaxBin/get_contig2bin.smk'

## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/contig2bin.tsv")
    include:
        '../../Binning/MetaBAT/metabat.smk'

## MetaDecoder
if "MetaDecoder" in BIN_STEPS:
    inputs_binning.append("Binning/MetaDecoder/contig2bin.tsv")
    include:
        '../../Binning/MetaDecoder/metadecoder.smk'

## SemiBin
if "SemiBin" in BIN_STEPS:
    inputs_binning.append("Binning/SemiBin/contig2bin.tsv")
    include:
        '../../Binning/SemiBin/semibin.smk'
## COMEBin
if "COMEBin" in BIN_STEPS:
    inputs_binning.append("Binning/COMEBin/contig2bin.tsv")
    include:
        '../../Binning/COMEBin/comebin.smk'

## GRiD
include:
    '../../Binning/GRiD/grid_dynamic.smk'

## master command
rule BINNING:
    input:
        ancient(inputs_binning)
    output:
        touch('status/binning.done')
