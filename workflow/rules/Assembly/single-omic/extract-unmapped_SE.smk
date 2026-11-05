SAMTOOLS_MEM = str(round(float(MAIN_BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['cluster_specs']['mem_unit']

EXTRACT_UNMAPPED_SHELL = """
                         echo "Mapping SE reads." >> {log} 2>&1
                         bwa-mem2 mem -v 1 -K {BWA_BATCH} \
                                      -t {threads} \
                                      {input[1]} \
                                      {input[0]} 2>> {log} | \
                         samtools fastq -0 {output[0]} \
                                        -f 4 \
                                        --threads 32 \
                                        -n \
                                        - 2>> {log}
                         """

rule extract_unmapped:
    input:
        'Preprocessing/{type}.se.preprocessed.fq.gz',
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.0123' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.amb' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.ann' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.bwt.2bit.64' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.pac' % IMP_ASSEMBLER,
    output:
        'Assembly/intermediary/{type}.se.unmapped.fq'
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_BIGMEMCORE_MB
    threads:
        getThreads(MAIN_BIGCORENO, bigmem=True)
    params:
        tmp_dir=TMPDIR
    conda:
        ENVDIR + "/IMP_mapping.yaml"
    log:
        "logs/assembly/extract_unmapped.{type}.log"
    message:
        "extract_unmapped: Extracting unmapped {wildcards.type} reads from megahit assembly."
    shadow: "copy-minimal"
    shell:
        EXTRACT_UNMAPPED_SHELL
