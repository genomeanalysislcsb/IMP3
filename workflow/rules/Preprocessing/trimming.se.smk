if config['nextseq']:
    rule trimming_se:
        input:
            'Preprocessing/{type}.se.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'Preprocessing/{type}.se.trimmoed.fq'
        params:
            runtime="12:00:00",
            mem = MEMCORE
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming_se.{type}.log"
        message: "trimming: Trimmming {wildcards.type} SE reads."
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

else:
    rule trimming_se:
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
            trimmomatic PE -threads {threads} {input[0]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """


