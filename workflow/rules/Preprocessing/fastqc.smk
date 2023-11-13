rule fastqc_raw:
    input:
        'Preprocessing/{type}.r1.fq',
        'Preprocessing/{type}.r2.fq'
    output:
        'Stats/{type}/{type}.r1_fastqc.html',
        'Stats/{type}/{type}.r2_fastqc.html',
        'Stats/{type}/{type}.r1_fastqc.zip',
        'Stats/{type}/{type}.r2_fastqc.zip',
    threads:
        getThreads(2)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        outputdir = 'Stats/{type}'
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log:
        "logs/preprocessing_fastqc_raw.{type}.log"
    message:
        "fastqc_raw: Running fastQC on raw {wildcards.type} reads."
    shell:
        """
        mkdir -p {params.outputdir}
        fastqc -o {params.outputdir} -f fastq {input[0]} {input[1]} -t {threads} -d {TMPDIR} > {log} 2>&1
        """

rule fastqc_preprocessed:
    input:
        'Preprocessing/{type}.r1.preprocessed.fq',
        'Preprocessing/{type}.r2.preprocessed.fq',
        'Preprocessing/{type}.se.preprocessed.fq',
        'Stats/{type}/{type}.r1_fastqc.html'
    output:
        'Stats/{type}/{type}.r1.preprocessed_fastqc.html',
        'Stats/{type}/{type}.r2.preprocessed_fastqc.html',
        'Stats/{type}/{type}.se.preprocessed_fastqc.html',
        'Stats/{type}/{type}.r1.preprocessed_fastqc.zip',
        'Stats/{type}/{type}.r2.preprocessed_fastqc.zip',
        'Stats/{type}/{type}.se.preprocessed_fastqc.zip'
    threads:
        getThreads(3)
    resources:
        runtime = VERY_SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        outputdir = 'Stats/{type}'
    log:
        "logs/preprocessing_fastqc_preprocessed.{type}.log"
    message:
        "fastqc_preprocessed: Running fastqc on preprocessed {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    shell:
        """
        if [[ -s {input[2]} ]]
        then
          fastqc -o {params.outputdir} -f fastq {input[0]} {input[1]} {input[2]} -t {threads} -d {TMPDIR} > {log} 2>&1
        else
          echo "{input[2]} is empty, FastQC will not run for that file" > {log} 2>&1
          fastqc -o {params.outputdir} -f fastq {input[0]} {input[1]} -t {threads} -d {TMPDIR} > {log} 2>&1
          echo "Creating empty files for {input[2]}..." > {log} 2>&1
          touch {output[2]}
        echo "No single end reads were generated" > {output[5]} > {log} 2>&1
        fi
        """
