# include rules for the Binning part of the workflow

include:
    'vizbin.smk'
include:
    'binny.smk'
include:
    'binny2dastool.smk'

## master command
rule BINNING:
    input:
        "Binning/binny/scaffold2bin.tsv"
    output:
        touch('status/binning.done')
