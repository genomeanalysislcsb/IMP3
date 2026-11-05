# include rules for the Analysis part of the workflow

include:
    '../../../Analysis/gene-depth.smk'
include:
    '../../../Analysis/bakta.smk'
include:
    '../../../Analysis/variant.smk'
include:
    '../../../Analysis/annotate.smk'
include:
    '../../../Proteomics/prot_db_prep.smk'

if len(MT) == 2:
    include:
        '../../../Analysis/featureCount.smk'
else:
    include:
        '../../../Analysis/featureCount_SE.smk'
    
# master command
rule ANALYSIS:
    input:
        expand("Stats/{type}/gene.thresholds.bed.gz", type=TYPES),
        expand("Stats/{type}/{type}.gene_depth.avg", type=TYPES),
        expand("Stats/{type}/{type}.gene.len", type=TYPES),
        expand("Analysis/snps/{type}.variants.samtools.vcf.gz", type=TYPES),
        expand("Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv",
               type=TYPES,
               db=["ID", "Name"],
               feature=["CRISPR", "CDS", "crispr-repeat", "crispr-spacer", "ncRNA", "oriC", "oriT", "regulatory_region", "rRNA", "tmRNA", "tRNA"]),
        aggregate_mantis_fc_per_contig
            if config["annotation"] == "mantis"
            else [],
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        expand("Analysis/annotation/{type}.{db}.{feature}_counts.tsv",
               type=TYPES,
               db=["ID", "Name"],
               feature=["CRISPR", "CDS", "crispr-repeat", "crispr-spacer", "ncRNA", "oriC", "oriT", "regulatory_region", "rRNA", "tmRNA", "tRNA"]),
        aggregate_mantis_fc
            if config["annotation"] == "mantis"
            else [],
        # "Analysis/annotation/intermediary.tar.gz" if 'binning' not in IMP_STEPS else [],
        "Analysis/annotation/proteomics.final.faa.gz"
    output:
        touch('status/analysis.done')
