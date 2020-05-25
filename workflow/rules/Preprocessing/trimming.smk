if config['nextseq']:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'Preprocessing/{type}.r1.trimmed.fq',
            'Preprocessing/{type}.se1.trimmed.fq',
            'Preprocessing/{type}.r2.trimmed.fq',
            'Preprocessing/{type}.se2.trimmed.fq'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """

    rule nextseq_trimming:
        input:
            'Preprocessing/{type}.{read}.trimmed.fq'
        output:
            'Preprocessing/{type}.{read}.trimmed.fq',
            'Preprocessing/{type}.{read}.trimmed.fq.gz'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_nextseq_trimming.{type}.{read}.log"
        message: "trimming: Trimmming poly-G's off {wildcards.type} {wildcards.read} reads."
        shell:
            """
            cutadapt --nextseq-trim={config[trimmomatic][trailing]} -o {output[0]} {input} >> {log} 2>&1
            gzip {input} >> {log} 2>&1
            """
else:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'Preprocessing/{type}.r1.trimmed.fq',
            'Preprocessing/{type}.se1.trimmed.fq',
            'Preprocessing/{type}.r2.trimmed.fq',
            'Preprocessing/{type}.se2.trimmed.fq'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """


rule cat_se_trimmed:
    input:
        'Preprocessing/{type}.se1.trimmed.fq',
        'Preprocessing/{type}.se2.trimmed.fq'
    output:
        'Preprocessing/{type}.se.trimmed.fq'
    threads: 1
    params:
        runtime="12:00:00",
        mem = MEMCORE
#    log: "logs/preprocessing_cat_se_trimmed.{type}.log"
    message: "cat_se_trimmed: Concatenating trimmed single end {wildcards.type} reads"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
#        gzip {input}
        """
