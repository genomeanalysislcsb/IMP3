# include rules for the Binning part of the workflow

# DAStool
include:
    '../../../Binning/DAStool/dastool.smk'

# MetaWrap
include:
    '../../../Binning/MetaWrap/metawrap.smk'

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']
inputs_binning = ["Binning/DASTool/selected_DASTool_summary.txt",
                  "Binning/MetaWrap/bins.stats"]

## MaxBin doesn't work without mg coverage
## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/scaffold2bin.tsv")
    include:
        '../../../Binning/MetaBAT/metabat.mt.smk'


## GRiD isn't run, but it contains a rule used by GTDBtk
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
