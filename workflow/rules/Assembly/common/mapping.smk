SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + "G"

rule mapping_on_assembly:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        "Assembly/%s.assembly.merged.fa.amb" % ASS,
        "Assembly/%s.assembly.merged.fa.bwt" % ASS,
        "Assembly/%s.assembly.merged.fa.pac" % ASS,
        "Assembly/%s.assembly.merged.fa.sa" % ASS,
        "Assembly/%s.assembly.merged.fa.ann" % ASS,
        'Assembly/%s.assembly.merged.fa' % ASS
    output:
        'Assembly/{type}.reads.sorted.bam'
    params:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: BIGCORENO
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_mapping_on_assembly.{type}.log"
    message: "mapping_on_assembly: Mapping {wildcards.type} reads on merged assembly."
    shell:
        """
        SAMHEADER="@RG\\tID:{SAMPLE}\\tSM:{wildcards.type}"
        PREFIX=Assembly/{wildcards.type}.reads
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

rule index_bam:
    input:
        '{sorted_bam}'
    output:
        '{sorted_bam}.bai'
    params:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_index.{sorted_bam}.log"
    message: "index_bam: Indexing {wildcards.sorted_bam}."
    shell:
        """
        samtools index {input} > {log} 2>&1
        """



