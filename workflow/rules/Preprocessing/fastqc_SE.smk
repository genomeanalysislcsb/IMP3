rule fastqc_raw:
    input:
        'Preprocessing/{type}.se.fq'
    output:
        'Stats/{type}/{type}.se_fastqc.html',
        'Stats/{type}/{type}.se_fastqc.zip'
    threads: getThreads(2)
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    params:
        outputdir = 'Stats/{type}'
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    log: "logs/preprocessing_fastqc_raw.{type}.log"
    message: "fastqc_raw: Running fastQC on raw {wildcards.type} reads."
    shell:
        """
        mkdir -p {params.outputdir}
        fastqc -o {params.outputdir} -f fastq {input[0]} -t {threads} -d {TMPDIR} > {log} 2>&1
        """

rule fastqc_preprocessed:
    input:
        'Preprocessing/{type}.se.preprocessed.fq',
        'Stats/{type}/{type}.se_fastqc.html'
    output:
        'Stats/{type}/{type}.se.preprocessed_fastqc.html',
        'Stats/{type}/{type}.se.preprocessed_fastqc.zip'
    threads: getThreads(3)
    resources:
        runtime = "1:00:00",
        mem = MEMCORE
    params:
        outputdir = 'Stats/{type}'
    log: "logs/preprocessing_fastqc_preprocessed.{type}.log"
    message: "fastqc_preprocessed: Running fastqc on preprocessed {wildcards.type} reads."
    conda:
        os.path.join(ENVDIR, "IMP_preprocessing.yaml")
    shell:
        """
        if [[ -s {input[0]} ]]
        then
          fastqc -o {params.outputdir} -f fastq {input[0]} -t {threads} -d {TMPDIR} > {log} 2>&1
        else
          echo "{input[0]} is empty, FastQC will not run for that file" > {log} 2>&1
          touch {output[0]}
          echo "No single end reads were generated" > {output[1]} > {log} 2>&1
        fi
        """
