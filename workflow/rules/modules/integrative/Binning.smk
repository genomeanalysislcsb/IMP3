# include rules for the Binning part of the workflow

# DAStoolBinning/DASTool/selected_DASTool
include:
    '../../Binning/DAStool/dastool.smk'

# MetaWrap
include:
    '../../Binning/MetaWrap/metawrap.smk'

mw_min_comp = config['binning']['MetaWrap']['min_completeness']
mw_max_cont = config['binning']['MetaWrap']['max_contamination']
# inputs_binning = ["status/binning_GRiD.done", "Binning/DASTool/selected_DASTool_summary.txt", "Binning/DASTool"
#                   "Binning/MetaWrap/bins.stats"]
inputs_binning = ["status/binning_GRiD.done", "Binning/MetaWrap/bins.stats", "Binning/DASTool"]

## abawaca
if "abawaca" in BIN_STEPS:
    inputs_binning.append("Binning/abawaca/contig2bin.tsv")
    inputs_binning.append("Binning/abawaca")
    include:
        '../../Binning/abawaca/abawaca.smk'

## CONCOCT
if "CONCOCT" in BIN_STEPS:
    inputs_binning.append("Binning/CONCOCT/contig2bin.tsv")
    inputs_binning.append("Binning/CONCOCT")
    include:
        '../../Binning/CONCOCT/concoct.smk'

## binny
if "binny" in BIN_STEPS:
    inputs_binning.append("Binning/binny/contig2bin.tsv")
    inputs_binning.append("Binning/binny")
    include:
        '../../Binning/binny/binny_subworkflow.smk'

## MaxBin
if "MaxBin" in BIN_STEPS:
    inputs_binning.append("Binning/MaxBin/maxbin_res.log")
    inputs_binning.append("Binning/MaxBin")
    include:
        '../../Binning/MaxBin/maxbin.smk'
    include:
        '../../Binning/MaxBin/get_contig2bin.smk'

## MetaBAT
if "MetaBAT" in BIN_STEPS:
    inputs_binning.append("Binning/MetaBAT/contig2bin.tsv")
    inputs_binning.append("Binning/MetaBAT")
    include:
        '../../Binning/MetaBAT/metabat.smk'

## MetaDecoder
if "MetaDecoder" in BIN_STEPS:
    inputs_binning.append("Binning/MetaDecoder/contig2bin.tsv")
    inputs_binning.append("Binning/MetaDecoder")
    include:
        '../../Binning/MetaDecoder/metadecoder.smk'

## SemiBin
if "SemiBin" in BIN_STEPS:
    inputs_binning.append("Binning/SemiBin/contig2bin.tsv")
    inputs_binning.append("Binning/SemiBin")
    include:
        '../../Binning/SemiBin/semibin.smk'
## COMEBin
if "COMEBin" in BIN_STEPS:
    inputs_binning.append("Binning/COMEBin/contig2bin.tsv")
    inputs_binning.append("Binning/COMEBin")
    include:
        '../../Binning/COMEBin/comebin.smk'

## GRiD
include:
    '../../Binning/GRiD/grid_dynamic.smk'

if not 'taxonomy' in IMP_STEPS:
    inputs_binning.append("Binning/per_bin_results.tar.gz")

## master command
rule BINNING:
    input:
        ancient(inputs_binning)
    output:
        touch('status/binning.done')
