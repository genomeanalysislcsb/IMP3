# include rules for the Binning part of the workflow

## DAStoolBinning/DASTool/selected_DASTool
include:
    '../../Binning/DAStool/dastool.smk'

# MetaWrap
include:
    '../../Binning/MetaWrap/metawrap.smk'

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']
inputs_binning = ["Binning/DASTool/selected_DASTool_summary.txt", "status/binning_GRiD.done",
                  f"Binning/MetaWrap/bins.stats"]

## abawaca
if "abawaca" in BIN_STEPS:
    inputs_binning.append("Binning/abawaca/scaffold2bin.tsv")
    include:
        '../../Binning/abawaca/abawaca.smk'

## CONCOCT
if "CONCOCT" in BIN_STEPS:
    inputs_binning.append("Binning/CONCOCT/scaffold2bin.tsv")
    include:
        '../../Binning/CONCOCT/concoct.smk'

## binny
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/scaffold2bin.tsv")
    include:
        '../../Binning/binny/binny_subworkflow.smk'

## MaxBin
if "MaxBin" in BIN_STEPS:
    inputs_binning.append("Binning/MaxBin/maxbin_res.log")
    include:
        '../../Binning/MaxBin/maxbin.smk'
    include:
        '../../Binning/MaxBin/get_contig2bin.smk'

## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/scaffold2bin.tsv")
    include:
        '../../Binning/MetaBAT/metabat.smk'

## MetaDecoder
if "MetaDecoder" in BIN_STEPS:
    inputs_binning.append("Binning/MetaDecoder/scaffold2bin.tsv")
    include:
        '../../Binning/MetaDecoder/metadecoder.smk'

## SemiBin
if "SemiBin" in BIN_STEPS:
    inputs_binning.append("Binning/SemiBin/scaffold2bin.tsv")
    include:
        '../../Binning/SemiBin/semibin.smk'

## GRiD
include:
    '../../Binning/GRiD/grid_dynamic.smk'

if not 'taxonomy' in IMP_STEPS:
    inputs_binning.append("Binning/per_bin_results.tar.gz")


## master command
rule BINNING:
    input:
        inputs_binning
    output:
        touch('status/binning.done')
