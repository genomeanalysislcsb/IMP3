# include rules for the Analysis part of the workflow

include:
    'prokkaC.smk'
include:
    'annotate.smk'

# master command
rule ANALYSIS:
    input:
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        "Analysis/annotation/intermediary.tar.gz",
    output:
        touch('status/analysis.done')
