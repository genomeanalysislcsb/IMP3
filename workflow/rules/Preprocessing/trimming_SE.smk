if config['nextseq']:
    rule trimming:
        input:
            'Preprocessing/{type}.se.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'Preprocessing/{type}.se.trimmoed.fq'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic SE -threads {threads} {input[0]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """

    rule nextseq_trimming:
        input:
            'Preprocessing/{type}.{read}.trimmoed.fq'
        output:
            'Preprocessing/{type}.{read}.trimmed.fq',
            'Preprocessing/{type}.{read}.trimmoed.fq.gz'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_nextseq_trimming.{type}.{read}.log"
        message: "trimming: Trimmming poly-G's off {wildcards.type} {wildcards.read} reads."
        shell:
            """
            cutadapt --nextseq-trim={config[trimmomatic][trailing]} -o {output} {input} >> {log} 2>&1
            gzip {input} >> {log} 2>&1
            """
else:
    rule trimming:
        input:
            'Preprocessing/{type}.se.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'Preprocessing/{type}.se.trimmed.fq'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic SE -threads {threads} {input[0]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """

