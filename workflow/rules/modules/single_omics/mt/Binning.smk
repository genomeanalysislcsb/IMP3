# include rules for the Binning part of the workflow

## DAStool
include:
    '../../../Binning/DAStool/dastool.smk'

inputs_binning = ["Binning/selected_DASTool_summary.txt"]

## MaxBin doesn't work without mg coverage
## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/metabat_res")
    include:
        '../../../Binning/MetaBAT/metabat.mt.smk'
## binny doesn't work without mg coverage
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/mt.vizbin.with-contig-names.points")
    include:
        '../../../Binning/binny/vizbin.smk'


## GRiD isn't run, but it contains a rule used by GTDBtk
include:
    '../../../Binning/GRiD/grid_dynamic.smk'

if not 'taxonomy' in IMP_STEPS:
    inputs_binning.append("Binning/per_bin_results.tar.gz")


## master command
rule BINNING:
    input:
        inputs_binning
    output:
        touch('status/binning.done')
