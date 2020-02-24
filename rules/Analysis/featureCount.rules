rule feature_count:
    input:
        "Assembly/{type}.reads.sorted.bam",
        "Analysis/annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    output:
        report("Analysis/annotation/{type}.{db}.counts.tsv",category="Analysis")
    params:
        lambda wildcards: config["featureCountsStranding"][wildcards.type],
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 8
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/analysis_feature_count.{type}.{db}.log"
    message: "feature_count: Quantifying {wildcards.db} with {wildcards.type} reads."
    shell:
        """
	featureCounts -p -O -t CDS -g {wildcards.db} -o {output} -s {params[0]} -a {input[1]} -T {threads} {input[0]} > {log} 2>&1
	"""

rule feature_count_gene:
    input:
        "Analysis/annotation/annotation.filt.gff",
        "Assembly/{type}.reads.sorted.bam"
    output:
        report("Analysis/annotation/{type}.{feature}_counts.tsv",category="Analysis")
    params:
        lambda wildcards: config["featureCountsStranding"][wildcards.type],
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 8
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/analysis_feature_count_gene.{feature}.{type}.genes.log"
    message: "feature_count_gene: Quantifying {wildcards.feature} with {wildcards.type} reads."
    shell:
        """
        featureCounts -p -O -t {wildcards.feature} -g ID -o {output} -s {params[0]}  -a {input[0]} -T {threads} {input[1]} > {log} 2>&1
        """
