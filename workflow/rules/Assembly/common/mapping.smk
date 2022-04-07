SAMTOOLS_MEM = str(round(float(BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['mem']['mem_unit']

rule mapping_on_assembly:
    input:
        r_1='Preprocessing/{type}.r1.preprocessed.fq',
        r_2='Preprocessing/{type}.r2.preprocessed.fq',
        r_se='Preprocessing/{type}.se.preprocessed.fq',
        assembly_amb="Assembly/%s.assembly.merged.fa.amb" % ASS,
        assembly_bwt="Assembly/%s.assembly.merged.fa.bwt" % ASS,
        assembly_pac="Assembly/%s.assembly.merged.fa.pac" % ASS,
        assembly_sa="Assembly/%s.assembly.merged.fa.sa" % ASS,
        assembly_ann="Assembly/%s.assembly.merged.fa.ann" % ASS,
        assembly='Assembly/%s.assembly.merged.fa' % ASS
    output:
        'Assembly/{type}.reads.sorted.bam'
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    threads: getThreads(BIGCORENO)
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_mapping_on_assembly.{type}.log"
    message: "mapping_on_assembly: Mapping {wildcards.type} reads on merged assembly."
    shell:
        """
        SAMHEADER="@RG\\tID:{SAMPLE}\\tSM:{wildcards.type}"
        PREFIX=Assembly/{wildcards.type}.reads
        # merge paired and se
        samtools merge --threads {threads} -f $PREFIX.merged.bam \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input.assembly} {input.r_1} {input.r_2} 2>> {log}| \
         samtools view --threads {threads} -bS -) \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input.assembly} {input.r_se} 2>> {log}| \
         samtools view --threads {threads} -bS -) 2>> {log}
        # sort
        samtools sort --threads {threads} -m {SAMTOOLS_MEM} $PREFIX.merged.bam > $PREFIX.sorted.bam 2>> {log}
        rm $PREFIX.merged.bam
        """
    # """
    # SAMHEADER="@RG\\tID:{SAMPLE}\\tSM:{wildcards.type}"
    # PREFIX=Assembly/{wildcards.type}.reads
    # # merge paired and se
    # samtools merge --threads {threads} -f $PREFIX.merged.bam \
    #  <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input.assembly} {input.r_1} {input.r_2} 2>> {log} | \
    #  samtools view --threads {threads} -q 10 -m 50 -F 4 -b -) \
    #  <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input.assembly} {input.r_se} 2>> {log} | \
    #  samtools view --threads {threads} -q 10 -m 50 -F 4 -b -) 2>> {log}
    # # sort
    # samtools sort --threads {threads} -m {resources.mem} $PREFIX.merged.bam | \
    #  samtools view --threads {threads} -f 256 -b - > $PREFIX.sorted.bam 2>> {log}
    # rm $PREFIX.merged.bam
    # """

rule index_bam:
    input:
        '{sorted_bam}'
    output:
        '{sorted_bam}.bai'
    resources:
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



