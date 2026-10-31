# Get raw read stats
rule get_raw_stats:
    input:
        raw = expand('Preprocessing/{{type}}.{read}.fq.gz', read=['se'])
    output:
        temp('Stats/{type}/{type}.raw.read_stats.tsv')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        12
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/analysis/{type}.raw.read_stats.log"
    message:
        "Getting stats for raw reads."
    shadow: "copy-minimal"
    shell:
        """
        echo "Getting stats for {input} ..." > {log}
        seqkit stats --threads {threads} \
                     --all \
                     --tabular \
                     --basename \
                     {input} > {output} 2>> {log}
        echo "Done." >> {log}
        """

# Count reads using seqit for each steps fastq files
rule count_step_reads:
    input:
        expand('Preprocessing/{{type}}.{read}.{{filterstep}}.fq.gz', read=['se'])
    output:
        temp('Stats/{type}/{type}.{filterstep}.read_stats.tsv')
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        12
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/analysis/{type}.{filterstep}.read_stats.log"
    message:
        "Getting stats for {wildcards.type} {wildcards.filterstep} reads."
    shadow: "minimal"
    shell:
        """
        echo "Counting reads in {input} ..." > {log}
        seqkit stats --threads {threads} \
                     --all \
                     --tabular \
                     --basename \
                     {input} > {output} 2>> {log}
        echo "Done." >> {log}
        """

rule gather_type_read_stats:
    input:
        type_read_stats = gather_filter_step_type_read_stats
    output:
        type_all_read_stats = report('Stats/{type}/{type}.read_stats.tsv',category="Preprocessing")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    log:
        "logs/analysis/{type}.read_stats.log"
    message:
        "reads_count: Getting length of {wildcards.type} fastq files."
    shadow: "copy-minimal"
    shell:
        """
        # Get header from first file
        echo "Getting header from first file ..." > {log}
        head -n 1 {input.type_read_stats[0]} > {output.type_all_read_stats} 2>> {log}

        # Add stats from all files, skippin headers
        for f in {input.type_read_stats}; do
            echo "Adding stats from ${{f}} ..." >> {log}
            tail -n +2 ${{f}} >> {output.type_all_read_stats} 2>> {log}
        done

        echo "Done." >> {log}
        """

rule reads_gzip:
    input:
        togz=gather_reads
    output:
        touch("status/preprocessing_{type}.reads.compressed")
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    message:
        "reads_gzip: Gzipped fastq files {input}."

# rule gzip_files:
#     input:
#         # Make sure its done at the very end
#         check = "status/workflow.done",
#         file = 'Preprocessing/{type}.{read}.preprocessed.fq.gz'
#     output:
#         "Preprocessing/{type}.{read}.preprocessed.fq.gz"
#     resources:
#         runtime = SHORT_RUNTIME,
#         mem_mb = MAIN_MEMCORE_MB
#     threads:
#         getThreads(16)
#     conda:
#         os.path.join(ENVDIR, "IMP_preprocessing.yaml")
#     log:
#         "logs/misc/compress.{type}.{read}.log" 
#     message:
#         "gzip_files: Gzipping fastq files {input}."
#     shadow: "minimal"
#     shell:
#         """
#         echo "Compressing {input.file} to {output}." > {log}
#         pigz -c -p {threads} \
#             -q \
#             {input.file} > {output} 2>> {log}
#         echo "Done." >> {log}
#        """


checkpoint reads_gziped:
    input:
        'Stats/{type}/{type}.read_stats.tsv'
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
