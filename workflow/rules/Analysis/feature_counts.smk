rule analysis_feature_counts:
    input:
        bam="Assembly/{omic}.reads.sorted.bam",
        gff="Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        "Analysis/annotation/{omic}.{db}.counts.tsv"
    log:
        "logs/analysis_feature_counts.{omic}.{db}.log"
    params:
        params=lambda wildcards: "-p -O -t CDS -g {} -s {}".format(wildcards.db, get_feature_counts_stranded(wildcards.omic))
    threads: 8
    resources:
        runtime="8:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "featureCounts: {input}"
    shell:
        "featureCounts -T {threads} {params.params} -o {output} -a {input.gff} {input.bam} &> {log} || touch {output} {output}.failed"

rule analysis_feature_counts_gene:
    input:
        bam="Assembly/{omic}.reads.sorted.bam",
        gff="Analysis/annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "Analysis/annotation/{omic}.{feature}_counts.tsv"
    log:
        "logs/analysis_feature_counts_gene.{omic}.{feature}.genes.log"
    wildcard_constraints:
        omic="|".join(get_omics()),
        # TODO: feature
    params:
        params=lambda wildcards: "-p -O -t {} -g ID -s {}".format(wildcards.feature, get_feature_counts_stranded(wildcards.omic))
    threads: 8
    resources:
        runtime="8:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_annotation.yaml")
    message:
        "featureCounts: {input}, {wildcards.feature}"
    shell:
        "featureCounts -T {threads} {params.params} -o {output} -a {input.gff} {input.bam} &> {log} || touch {output} {output}.failed"
