rule fastqc_raw:
    input:
        reads_r1_raw = 'Preprocessing/{type}.r1.fq.gz',
        reads_r2_raw = 'Preprocessing/{type}.r2.fq.gz'
    output:
        'Stats/{type}/{type}.r1_fastqc.html',
        'Stats/{type}/{type}.r2_fastqc.html',
        'Stats/{type}/{type}.r1_fastqc.zip',
        'Stats/{type}/{type}.r2_fastqc.zip',
    threads:
        getThreads(8)
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
    shadow: "minimal"
    shell:
        """
        mkdir -p {params.out_dir} > {log} 2>&1
        fastqc --outdir {params.out_dir} \
               --format fastq \
               --threads {threads} \
               --dir {params.tmp_dir} \
               {input.reads_r1_raw} \
               {input.reads_r2_raw} >> {log} 2>&1
        """

rule fastqc_preprocessed:
    input:
        reads_r1 = 'Preprocessing/{type}.r1.preprocessed.fq.gz',
        reads_r2 = 'Preprocessing/{type}.r2.preprocessed.fq.gz',
        reads_se = 'Preprocessing/{type}.se.preprocessed.fq.gz',
        fqc_r1_raw_report = 'Stats/{type}/{type}.r1_fastqc.html'
    output:
        fqc_r1_report = 'Stats/{type}/{type}.r1.preprocessed_fastqc.html',
        fqc_r2_report = 'Stats/{type}/{type}.r2.preprocessed_fastqc.html',
        fqc_se_report = 'Stats/{type}/{type}.se.preprocessed_fastqc.html',
        fqc_r1_output = 'Stats/{type}/{type}.r1.preprocessed_fastqc.zip',
        fqc_r2_output = 'Stats/{type}/{type}.r2.preprocessed_fastqc.zip',
        fqc_se_output = 'Stats/{type}/{type}.se.preprocessed_fastqc.zip'
    threads:
        getThreads(8)
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
    shadow: "minimal"
    shell:
        """
        if [[ -s {input.reads_se} ]]; then
          fastqc --outdir {params.out_dir} \
                 --format fastq \
                 --threads {threads} \
                 --memory {resources.mem_mb} \
                 --dir {params.tmp_dir} \
                 {input.reads_r1} \
                 {input.reads_r2} \
                 {input.reads_se} > {log} 2>&1
        else
          echo "{input.reads_se} is empty, FastQC will not run for that file" > {log} 2>&1
          fastqc --outdir {params.out_dir} \
                 --format fastq \
                 --threads {threads} \
                 --memory {resources.mem_mb} \
                 --dir {params.tmp_dir} \
                 {input.reads_r1} \
                 {input.reads_r2} >> {log} 2>&1
          echo "Creating empty files for {input.reads_se}..." >> {log} 2>&1
          touch {output.fqc_se_report}
          # TODO: This message is written to a zip file ???
          echo "No single end reads were generated" > {output.fqc_se_output} 2>> {log}
        fi
        """
