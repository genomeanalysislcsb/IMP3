SAMTOOLS_MEM = str(round(float(MAIN_BIGMEMCORE[:-1]) * 0.75 - 0.5)) + config['cluster_specs']['mem_unit']

EXTRACT_UNMAPPED_SHELL = """
                         export TMPDIR="{params.tmp_dir}"
                         TMP_FILE=$(mktemp "alignment_XXXXXX.bam")

                         echo "Mapping PE reads." >> {log} 2>&1

                         bwa-mem2 mem -v 1 -K {BWA_BATCH} \
                                      -t {threads} \
                                      {input[3]} \
                                      {input[0]} \
                                      {input[1]} 2>> {log} | \
                           samtools view --threads {threads} \
                                         -bS \
                                         - > ${{TMP_FILE}} 2> {log}
                                         
                         samtools merge --threads {threads} \
                                        -u \
                                        - \
                                        <(samtools view --threads {threads} \
                                                        -u \
                                                        -f 4 \
                                                        -F 264 \
                                                        ${{TMP_FILE}} 2>> {log}) \
                                        <(samtools view --threads {threads} \
                                                        -u \
                                                        -f 8 \
                                                        -F 260 \
                                                        ${{TMP_FILE}} 2>> {log}) \
                                        <(samtools view --threads {threads} \
                                                        -u \
                                                        -f 12 \
                                                        -F 256 \
                                                        ${{TMP_FILE}} 2>> {log}) 2>> {log} | \
                           samtools view --threads {threads} \
                                         -bF 0x800 \
                                         - 2>> {log} | \
                           samtools sort --threads {threads} \
                                         -m {SAMTOOLS_MEM} \
                                         -n \
                                         - 2>> {log} | \
                           samtools fastq -1 {output[0]} \
                                          -2 {output[1]} \
                                          -0 /dev/null \
                                          -s /dev/null \
                                          -f 4 \
                                          --threads {threads} \
                                          -n \
                                          - >> {log} 2>&1
                         
                         if [[ -s {input[2]} ]]; then
                           echo "Mapping SE reads." >> {log} 2>&1
                           bwa-mem2 mem -v 1 -K {BWA_BATCH} \
                                        -t {threads} \
                                        {input[3]} \
                                        {input[2]} 2>> {log} | \
                           samtools fastq -0 {output[2]} \
                                          -f 4 \
                                          --threads 32 \
                                          -n \
                                          - 2>> {log}
                         else
                           echo "There are no singletons reads. {input[2]} is empty. Generating empty {output[2]}" >> {log} 
                           touch {output[2]}
                         fi

                         rm -rf ${{TMP_FILE}}
                         """

rule extract_unmapped:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq.gz',
        'Preprocessing/{type}.r2.preprocessed.fq.gz',
        'Preprocessing/{type}.se.preprocessed.fq.gz',
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.0123' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.amb' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.ann' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.bwt.2bit.64' % IMP_ASSEMBLER,
        'Assembly/intermediary/{type}.%s_preprocessed.1.fa.pac' % IMP_ASSEMBLER,
    output:
        'Assembly/intermediary/{type}.r1.unmapped.fq',
        'Assembly/intermediary/{type}.r2.unmapped.fq',
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
