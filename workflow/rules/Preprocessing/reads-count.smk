rule reads_count:
    input:
        reads=reads_input_files
    output:
        report('Stats/{type}/{type}.read_counts.txt',category="Preprocessing")
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
#    log: "logs/analysis_read_numbers.{type}.log"
    message: "reads_count: Getting length of {wildcards.type} fastq files."
    run:
        for fs in {input.reads}:
            for f in fs:
                shell('echo -e "{f}\t$(($(wc -l < {f})/4))" >> {output}')

rule reads_gzip:
    input:
        togz=gather_reads
    output:
        "status/preprocessing_{type}.reads.compressed"
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    message: "reads_gzip: Gzipped fastq files {input}."
    shell:
        """
        touch {output}
        """

rule gzip_ori_files:
    input:
        'Preprocessing/{type}.{read}.fq'
    output:
        'Preprocessing/{type}.{read}.fq.gz'
    wildcard_constraints:
        read="r.",
        type="m."
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 5
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message: "gzip_ori_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} -q {input} || touch {input}"


rule gzip_files:
    input:
        'Preprocessing/{type}.{read}.{stage}.fq'
    output:
        "Preprocessing/{type}.{read}.{stage}.fq.gz"
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 5
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    message: "gzip_files: Gzipping fastq files {input}."
    shell:
        "pigz -p {threads} -q {input}"

checkpoint reads_gziped:
    input:
        'Stats/{type}/{type}.read_counts.txt'
    output:
        "status/preprocessing_{type}.reads.detected"
    params:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    message: "reads_gziped: Collecting reads for compression."
    run:
        shell("touch {output}")
        
