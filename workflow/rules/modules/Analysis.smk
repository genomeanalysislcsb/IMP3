include:
    "../Analysis/annotate.smk"
include:
    "../Analysis/feature_counts.smk"
include:
    "../Analysis/gene_depth.smk"
include:
    "../Analysis/misc.smk"
include:
    "../Analysis/prep_proteomics.smk"

rule ANALYSIS_READY:
    input:
        # GFF
        "Analysis/annotation/annotation_CDS_RNA_hmms.gff",
        # gene stats
        expand(
            "Stats/{omic}/annotation/{omic}.{ext}",
            omic=get_omics(), ext=["gene_depth.hist", "gene_depth.avg", "gene.len"]
        ),
        # variants
        expand(
            "Analysis/snps/{omic}.variants.samtools.vcf.gz",
            omic=get_omics(),
        ),
        # feature counts
        expand(
            "Analysis/annotation/{omic}.{feature}_counts.tsv",
            omic=get_omics(), feature=["rRNA","CDS"]
        ) if "mg" in get_omics() else [
            "Analysis/annotation/mt.CDS_counts.tsv"
        ],
        # HMM counts
        expand(
            "Analysis/annotation/{omic}.{db}.counts.tsv",
            omic=get_omics(), db=get_hmms()
        ),
        # proteomics
        "Analysis/annotation/proteomics.final.faa"
    output:
        touch("status/analysis.ready")

rule ANALYSIS:
    input:
        "status/analysis.ready",
        "Analysis/annotation/intermediary.tar.gz",
    output:
        touch("status/analysis.done")
