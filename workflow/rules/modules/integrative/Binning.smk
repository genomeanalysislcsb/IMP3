# include rules for the Binning part of the workflow

## DAStool
include:
    '../../Binning/DAStool/dastool.smk'

inputs_binning = ["Binning/selected_DASTool_summary.txt","status/binning_GRiD.done"]

## MaxBin
if "MaxBin" in BIN_STEPS:
    inputs_binning.append("Binning/MaxBin/maxbin_res.log")
    include:
        '../../Binning/MaxBin/maxbin.smk'
    include:
        '../../Binning/MaxBin/maxbin_get_depth.smk'
    include:
        '../../Binning/MaxBin/get_contig2bin.smk'
## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/metabat_res")
    include:
        '../../Binning/MetaBAT/metabat.smk'
## binny
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/scaffold2bin.tsv")
    include:
        '../../Binning/binny/vizbin.smk'
    include:
        '../../Binning/binny/binny.smk'
    include:
        '../../Binning/binny/binny2dastool.smk'
## GRiD
include:
    '../../Binning/GRiD/grid_dynamic.smk'

if not 'taxonomy' in IMP_STEPS:
    inputs_binning.append("Binning/selected_DASTool_bins/per_bin_results.tar.gz")


## master command
rule BINNING:
    input:
        inputs_binning
    output:
        touch('status/binning.done')
