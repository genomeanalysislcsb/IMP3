rule reads_count:
    input:
        reads=reads_input_files_SE
    output:
        report('Stats/{type}/{type}.read_counts.txt',category="Preprocessing")
    resources:
        runtime = "02:00:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    # log:
    #     "logs/analysis_read_numbers.{type}.log"
    message:
        "reads_count: Getting length of {wildcards.type} fastq files."
    run:
        for fs in {input.reads}:
            for f in fs:
                shell('echo -e "{f}\t$(($(wc -l < {f})/4))" >> {output}')

rule reads_gzip:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        togz=gather_reads
    output:
        "status/preprocessing_{type}.reads.compressed"
    resources:
        runtime = "00:01:00",
        mem_mb = MEMCORE_MB
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
        "status/workflow.done",
        'Preprocessing/{type}.{read}.fq'
    output:
        'Preprocessing/{type}.{read}.fq.gz'
    wildcard_constraints:
        read="se",
        type="mg"
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    threads:
        getThreads(5)
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "gzip_ori_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} \
              -q \
              {input} || touch {input}"


rule gzip_files:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        'Preprocessing/{type}.{read}.{stage}.fq'
    output:
        "Preprocessing/{type}.{read}.{stage}.fq.gz"
    resources:
        runtime = "04:00:00",
        mem_mb = MEMCORE_MB
    threads:
        getThreads(5)
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message:
        "gzip_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} \
              -q \
              {input}"

checkpoint reads_gziped:
    input:
        # Make sure its done at the very end
        "status/workflow.done",
        'Stats/{type}/{type}.read_counts.txt'
    output:
        "status/preprocessing_{type}.reads.detected"
    resources:
        runtime = "00:01:00",
        mem_mb = MEMCORE_MB
    threads:
        1
    message:
        "reads_gziped: Collecting reads for compression."
    run:
        shell("touch {output}")
        
