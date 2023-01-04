if config['nextseq']:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq'
        output:
            'Preprocessing/{type}.r1.trimmoed.fq',
            'Preprocessing/{type}.se1.trimmoed.fq',
            'Preprocessing/{type}.r2.trimmoed.fq',
            'Preprocessing/{type}.se2.trimmoed.fq'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = MEMCORE,
        params:
            adapter = lambda wildcards: config["trimmomatic"]["adapter"][wildcards.type]
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{params.adapter}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """

    rule nextseq_trimming_se:
        input:
            'Preprocessing/{type}.{read}.trimmoed.fq'
        output:
            'Preprocessing/{type}.{read}.trimmed.fq',
            'Preprocessing/{type}.{read}.trimmoed.fq.gz'
        wildcard_constraints:
            read="se[1|2]"
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = MEMCORE
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log: "logs/preprocessing_nextseq_trimming_se.{type}.{read}.log"
        message: "trimming: Trimmming poly-G's off {wildcards.type} {wildcards.read} reads."
        shell:
            """
            cutadapt -j {threads} --nextseq-trim={config[trimmomatic][trailing]} -o {output[0]} \
             --minimum-length {config[trimmomatic][minlen]} {input} >> {log} 2>&1
            pigz -p {threads} {input} >> {log} 2>&1
            """

    rule nextseq_trimming_pe:
        input:
            'Preprocessing/{type}.r1.trimmoed.fq',
            'Preprocessing/{type}.r2.trimmoed.fq'
        output:
            'Preprocessing/{type}.r1.trimmed.fq',
            'Preprocessing/{type}.r2.trimmed.fq',
            'Preprocessing/{type}.r1Cut.trimmed.fq',
            'Preprocessing/{type}.r2Cut.trimmed.fq',
            'Preprocessing/{type}.se1Cut.trimmed.fq',
            'Preprocessing/{type}.se2Cut.trimmed.fq',
            'Preprocessing/{type}.r1.trimmoed.fq.gz',
            'Preprocessing/{type}.r2.trimmoed.fq.gz'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = MEMCORE
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log: "logs/preprocessing_nextseq_trimming_pe.{type}.pe.log"
        message: "trimming: Trimmming poly-G's off {wildcards.type} pe reads."
        shell:
            """
            cutadapt -j {threads} --nextseq-trim={config[trimmomatic][trailing]} \
             --minimum-length {config[trimmomatic][minlen]}:{config[trimmomatic][minlen]} \
             -o {output[0]} -p {output[1]} --too-short-output {output[2]} --too-short-paired-output {output[3]} \
             {input[0]} {input[1]} >> {log} 2>&1
            cutadapt -j {threads} \
             --minimum-length {config[trimmomatic][minlen]} \
             -o {output[4]} {output[2]} >> {log} 2>&1
            cutadapt -j {threads} \
             --minimum-length {config[trimmomatic][minlen]} \
             -o {output[5]} {output[3]} >> {log} 2>&1
            pigz -p {threads} {input} >> {log} 2>&1
            """
    rule cat_se_trimmed:
        input:
            'Preprocessing/{type}.se1.trimmed.fq',
            'Preprocessing/{type}.se2.trimmed.fq',
            'Preprocessing/{type}.se1Cut.trimmed.fq',
            'Preprocessing/{type}.se2Cut.trimmed.fq',
            'Preprocessing/{type}.r1Cut.trimmed.fq',
            'Preprocessing/{type}.r2Cut.trimmed.fq'
        output:
            'Preprocessing/{type}.se.trimmed.fq',
        threads: getThreads(2)
        resources:
            runtime="12:00:00",
            mem = MEMCORE
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message: "cat_se_trimmed: Concatenating trimmed single end {wildcards.type} reads"
        shell:
            """
            cat {input[0]} {input[1]} {input[2]} {input[3]} > {output}
            pigz -p {threads} {input}
            """
else:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq'
        output:
            'Preprocessing/{type}.r1.trimmed.fq',
            'Preprocessing/{type}.se1.trimmed.fq',
            'Preprocessing/{type}.r2.trimmed.fq',
            'Preprocessing/{type}.se2.trimmed.fq'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = MEMCORE,
        params:
            adapter = lambda wildcards: config["trimmomatic"]["adapter"][wildcards.type]
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log: "logs/preprocessing_trimming.{type}.log"
        message: "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{params.adapter}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
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
        resources:
            runtime="12:00:00",
            mem = MEMCORE
        message: "cat_se_trimmed: Concatenating trimmed single end {wildcards.type} reads"
        shell:
            """
            cat {input[0]} {input[1]} > {output}
            """
