rule mapping_on_catalogue:
    input:
        lambda wildcards: expand("{raw_directory}/Preprocessing/{{type}}.{read}.preprocessed.fq",
                raw_directory=samples.loc[samples["sample"] == wildcards.sample,"path"],
                read=["r1","r2","se"]),
        expand("Catalogue/catalogue.{ext}",ext=["fasta.amb","fasta.bwt","fasta.pac","fasta.sa","fasta.ann","fasta"])
    output:
        'Catalogue/{sample}.{type}.reads.sorted.bam'
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/catalogue_mapping_on_catalogue.{type}.{sample}.log"
    message: "mapping_on_catalogue: Mapping {wildcards.type} reads on catalogue."
    shell:
        """
        SAMHEADER="@RG\\tID:{wildcards.sample}\\tSM:MG"
        PREFIX=Catalogue/{wildcards.sample}.{wildcards.type}.reads
        # merge paired and se
        samtools merge --threads {threads} -f $PREFIX.merged.bam \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[0]} {input[1]} 2>> {log}| \
         samtools view --threads {threads} -bS -) \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[2]} 2> {log}| \
         samtools view --threads {threads} -bS -) 2>> {log}
        # sort
        samtools sort --threads {threads} -m {SAMTOOLS_MEM} $PREFIX.merged.bam > $PREFIX.sorted.bam 2>> {log}
        rm $PREFIX.merged.bam
        """

rule feature_count_cat:
    input:
        gff = "Catalogue/catalogue.gff",
        bam = lambda wildcards: expand('Catalogue/{sample}.{{type}}.reads.sorted.bam',
              sample=samples.loc[samples[wildcards.type] == 1,"sample"])
    output:
        "Catalogue/{type}.counts.tsv"
    params:
        lambda wildcards: config["catalogue"]["mapping"]["featureCountsStranding"][wildcards.type]
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: getThreads(8)
    conda: ENVDIR + "/IMP_annotation.yaml"
    log: "logs/catalogue_feature_count_cat.{type}.log"
    message: "feature_count_cat: Quantifying {wildcards.type} reads on the catalogue."
    shell:
        """
        featureCounts -M --fraction -p --largestOverlap -t gene -g ID -o {output} -s {params[0]} -a {input.gff} -T {threads} {input.bam} > {log} 2>&1
        """
