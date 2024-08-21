#########################
### hmmer db counting ###
#########################

# TODO: Find a way to dynamically generate log names, since params are not accessible by log

# paired total feature sum
rule feature_count_hmmer:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        "Analysis/annotation/{type}.hmmer.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "total",
        count_type = "hmmer"
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(12)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis/featurecount/hmmer/total.{type}.{db}.paired.log"
    message:
        "feature_count: Quantifying {params.db} with {wildcards.type} {params.read_type} reads ({params.feature_sum} counts)."
    shell:
        """
        featureCounts {params.flags} \
                      -t {params.feature} \
                      -g {params.db} \
                      -o {output} \
                      -s {params.stranding} \
                      -a {input.gff} \
                      -T {threads} \
                      {input.map}> {log} 2>&1 \
                      || touch {output} \
                               status/analysis_featureCounts.{wildcards.type}.{params.db}.{params.read_type}.impossible
        """


# single per contig sum
use rule feature_count_hmmer as feature_count_hmmer_per_contig with:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        "Analysis/annotation/per_contig_counts/{type}.hmmer.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "per_contig",
        count_type = "hmmer"
    log:
        "logs/analysis/featurecount/hmmer/per_contig.{type}.{db}.single.log"

##########################
### Mantis db counting ###
##########################

# single total feature sum
use rule feature_count_hmmer as feature_count_mantis with:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.{db}.gff"
    output:
        "Analysis/annotation/{type}.mantis.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "total"
    log:
        "logs/analysis/featurecount/mantis/total.{type}.{db}.single.log"

# single per contig sum
use rule feature_count_hmmer as feature_count_mantis_per_contig with:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.{db}.gff"
    output:
        "Analysis/annotation/per_contig_counts/{type}.mantis.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = "CDS",
        read_type = "single",
        feature_sum = "per_contig"
    log:
        "logs/analysis/featurecount/mantis/per_contig.{type}.{db}.single.log"

##########################
### Gene level counting ##
##########################

# single total feature sum
use rule feature_count_hmmer as feature_count_gene with:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/{type}.{db}.{feature}_counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "single",
        feature_sum = "total",
        count_type = "gene"
    log:
        "logs/analysis/featurecount/gene/total.{type}.{db}.{feature}.single.log"

# single per contig sum
use rule feature_count_hmmer as feature_count_gene_per_contig with:
    input:
        map = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/per_contig_counts/{type}.{db}.{feature}_per_contig.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        db = lambda wildcards: wildcards.db,
        flags = "-f -O --fraction -M --fraction",
        feature = lambda wildcards: wildcards.feature,
        read_type = "single",
        feature_sum = "per_contig",
        count_type = "gene"
    log:
        "logs/analysis/featurecount/gene/per_contig.{type}.{db}.{feature}.single.log"
