# include rules for the Analysis part of the workflow

include:
    '../../Analysis/gene-depth.smk'
include:
    '../../Analysis/prokkaC.smk'
include:
    '../../Analysis/variant.smk'
include:
    '../../Analysis/annotate.smk'
include:
    '../../Analysis/prepProteomics.smk'
include:
    '../../Analysis/featureCount.smk'

# master command
rule ANALYSIS:
    input:
        expand("Stats/{type}/annotation/{type}.gene_depth.hist", type=TYPES),
        expand("Stats/{type}/annotation/{type}.gene_depth.avg", type=TYPES),
        expand("Stats/{type}/annotation/{type}.gene.len", type=TYPES),
        expand("Analysis/snps/{type}.variants.samtools.vcf.gz", type=TYPES),
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        expand("Analysis/annotation/{type}.{feature}_counts.tsv",type=TYPES,feature=["rRNA","CDS"]),
        expand("Analysis/annotation/{type}.{db}.counts.tsv",type=TYPES,db=config["hmm_DBs"].split()),
        "Analysis/annotation/intermediary.tar.gz",
        "Analysis/annotation/proteomics.final.faa"
    output:
        touch('status/analysis.done')
