rule reads_count:
    input:
        reads=reads_input_files_SE
    output:
        report('Stats/{type}/{type}.read_stats.tsv',category="Preprocessing")
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
        togz = gather_reads
    output:
        touch("status/preprocessing_{type}.reads.compressed")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "reads_gzip: Gzipped fastq files {input}."

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
    log:
        "logs/misc/compress.{type}.{read}.log" 
    message:
        "gzip_files: Gzipping fastq files {input}."
    shell:
        """
        echo "Compressing {input.file} to {output}." > {log}
        pigz -p {threads} \
             -q \
             {input.file} >> {log} 2>&1
        echo "Done." >> {log}
        """

checkpoint reads_gziped:
    input:
        file = 'Stats/{type}/{type}.read_stats.tsv'
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
        
