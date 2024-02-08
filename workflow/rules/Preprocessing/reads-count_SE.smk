rule reads_count:
    input:
        reads=reads_input_files_SE
    output:
        report('Stats/{type}/{type}.read_counts.txt',category="Preprocessing")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    # log:
    #     "logs/analysis/read_numbers.{type}.log"
    message:
        "reads_count: Getting length of {wildcards.type} fastq files."
    run:
        for fs in {input.reads}:
            for f in fs:
                shell('echo -e "{f}\t$(($(wc -l < {f})/4))" >> {output}')

rule reads_gzip:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        togz = gather_reads
    output:
        "status/preprocessing_{type}.reads.compressed"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "reads_gzip: Gzipped fastq files {input}."
    shell:
        """
        touch {output}
        """

rule gzip_ori_files:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        file = 'Preprocessing/{type}.{read}.fq'
    output:
        'Preprocessing/{type}.{read}.fq.gz'
    wildcard_constraints:
        read="se",
        type="mg"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(5)
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "gzip_ori_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} \
              -q \
              {input.file} || touch {input.file}"


rule gzip_files:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        file = 'Preprocessing/{type}.{read}.preprocessed.fq'
    output:
        "Preprocessing/{type}.{read}.preprocessed.fq.gz"
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(5)
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "gzip_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} \
              -q \
              {input.file}"

checkpoint reads_gziped:
    input:
        # Make sure its done at the very end
        check = "status/workflow.done",
        file = 'Stats/{type}/{type}.read_counts.txt'
    output:
        "status/preprocessing_{type}.reads.detected"
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "reads_gziped: Collecting reads for compression."
    run:
        shell("touch {output}")
        
