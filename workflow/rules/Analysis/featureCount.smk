rule feature_count_hmmer:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        "Analysis/annotation/{type}.hmmer.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O'
    resources:
        runtime = "04:00:00",
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(8)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis_feature_count.{type}.{db}.log"
    message:
        "feature_count: Quantifying {wildcards.db} with {wildcards.type} reads."
    shell:
        """
        featureCounts {params.flags} \
                      -t CDS \
                      -g {wildcards.db} \
                      -o {output} \
                      -s {params.stranding} \
                      -a {input.gff} \
                      -T {threads} \
                      {input.bam} > {log} 2>&1 \
                      || touch {output} \
                               status/analysis_featureCounts.{wildcards.type}.{wildcards.db}.impossible
        """

rule feature_count_mantis:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/intermediary/mantis/consensus_annotation.mantis.{db}.gff"
    output:
        "Analysis/annotation/{type}.mantis.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O'
    resources:
        runtime = "04:00:00",
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(8)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis_feature_count.{type}.mantis.{db}.log"
    message:
        "feature_count: Quantifying mantis {wildcards.db} with {wildcards.type} reads."
    shell:
        """
        featureCounts {params.flags} \
                      -t CDS \
                      -g {wildcards.db} \
                      -o {output} \
                      -s {params.stranding} \
                      -a {input.gff} \
                      -T {threads} \
                      {input.bam} > {log} 2>&1 \
                      || touch {output} \
                               status/analysis_featureCounts.{wildcards.type}.{wildcards.db}.impossible
            """


rule feature_count_gene:
    input:
        bam = "Assembly/{type}.reads.sorted.bam",
        gff = "Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/{type}.{feature}_counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O'
    resources:
        runtime = "04:00:00",
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(8)
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    log:
        "logs/analysis_feature_count_gene.{feature}.{type}.genes.log"
    message:
        "feature_count_gene: Quantifying feature {wildcards.feature} with {wildcards.type} reads."
    shell:
        """
        featureCounts {params.flags} \
                      -t {wildcards.feature} \
                      -g ID \
                      -o {output} \
                      -s {params.stranding} \
                      -a {input.gff} \
                      -T {threads} \
                      {input.bam} > {log} 2>&1 \
                      || touch {output} \
                               status/analysis_featureCounts.{wildcards.type}.{wildcards.feature}.impossible
        """

use rule feature_count_hmmer as feature_count_hmmer_per_contig with:
    output:
        "Stats/{type}/annotation/{type}.hmmer.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O -f'
    log:
        "logs/stats_feature_count_per_contig.{type}.{db}.log"
    message:
        "feature_count: Quantifying {wildcards.db} per contig with {wildcards.type} reads."


use rule feature_count_mantis as feature_count_mantis_per_contig with:
    output:
        "Stats/{type}/annotation/{type}.mantis.{db}.counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O -f'
    log:
        "logs/stats_feature_count_per_contig.{type}.mantis.{db}.log"
    message:
        "feature_count: Quantifying mantis {wildcards.db} per contig with {wildcards.type} reads."


use rule feature_count_gene as feature_count_gene_per_contig with:
    output:
        "Stats/{type}/annotation/{type}.{feature}_counts.tsv"
    params:
        stranding = lambda wildcards: config["featureCountsStranding"][wildcards.type],
        flags = '-p -O -f'
    log:
        "logs/stats_feature_count_gene_per_contig.{feature}.{type}.genes.log"
    message:
        "feature_count_gene: Quantifying feature {wildcards.feature} per contig with {wildcards.type} reads."