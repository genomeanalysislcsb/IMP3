# include rules for the Analysis part of the workflow

include:
    '../../../Analysis/gene-depth.smk'
include:
    '../../../Analysis/prokkaC.smk'
include:
    '../../../Analysis/variant.smk'
include:
    '../../../Analysis/annotate.smk'
include:
    '../../../Analysis/prepProteomics.smk'
include:
    '../../../Analysis/featureCount.smk'

# master command
rule ANALYSIS:
    input:
        expand("Stats/{type}/annotation/{type}.gene_depth.hist", type=TYPES),
        expand("Stats/{type}/annotation/{type}.gene_depth.avg", type=TYPES),
        expand("Stats/{type}/annotation/{type}.gene.len", type=TYPES),
        expand("Analysis/snps/{type}.variants.samtools.vcf.gz", type=TYPES),
        expand("Stats/{type}/annotation/{type}.hmmer.{db}.counts.tsv",type=TYPES,db=config["hmm_DBs"].split()) if config["annotation"] == "hmmer" else [],
        aggregate_mantis_fc_per_contig if config["annotation"] == "mantis" else [],
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        "Analysis/annotation/mt.CDS_counts.tsv",
        expand("Analysis/annotation/{type}.hmmer.{db}.counts.tsv",type=TYPES,db=config["hmm_DBs"].split()) if config["annotation"] == "hmmer" else [],
        aggregate_mantis_fc if config["annotation"] == "mantis" else [],
        # "Analysis/annotation/intermediary.tar.gz" if 'binning' not in IMP_STEPS else [],
        "Analysis/annotation/proteomics.final.faa"
    output:
        touch('status/analysis.done')
