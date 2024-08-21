rule fastqc_raw:
    input:
        reads_se_raw = 'Preprocessing/{type}.se.fq'
    output:
        'Stats/{type}/{type}.se_fastqc.html',
        'Stats/{type}/{type}.se_fastqc.zip'
    threads:
        getThreads(4)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        out_dir = 'Stats/{type}',
        tmp_dir = TMPDIR
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/preprocessing/fastqc_raw.{type}.log"
    message:
        "fastqc_raw: Running fastQC on raw {wildcards.type} reads."
    shell:
        """
        mkdir -p {params.out_dir} > {log} 2>&1
        fastqc --outdir {params.out_dir} \
               --format fastq \
               --threads {threads} \
               --memory {resources.mem_mb} \
               --dir {params.tmp_dir} \
               {input.reads_se_raw} >> {log} 2>&1
        """

rule fastqc_preprocessed:
    input:
        reads_se = 'Preprocessing/{type}.se.preprocessed.fq',
        fqc_raw_report = 'Stats/{type}/{type}.se_fastqc.html'
    output:
        fqc_se_report = 'Stats/{type}/{type}.se.preprocessed_fastqc.html',
        fqc_se_output = 'Stats/{type}/{type}.se.preprocessed_fastqc.zip'
    threads:
        getThreads(4)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        out_dir = 'Stats/{type}',
        tmp_dir = TMPDIR
    log:
        "logs/preprocessing/fastqc_preprocessed.{type}.log"
    message:
        "fastqc_preprocessed: Running fastqc on preprocessed {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    shell:
        """
        if [[ -s {input.reads_se} ]]; then
          fastqc --outdir {params.out_dir} \
                 --format fastq \
                 --threads {threads} \
                 --memory {resources.mem_mb} \
                 --dir {params.tmp_dir} \
                 {input.reads_se} > {log} 2>&1
        else
          echo "{input.reads_se} is empty, FastQC will not run for that file." >> {log} 2>&1
          touch {output.fqc_se_report}
          # TODO: This message is written to a zip file ???
          echo "No single end reads were generated" > {output.fqc_se_output} 2>> {log}
        fi
        """
