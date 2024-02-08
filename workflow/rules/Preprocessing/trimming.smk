import math

if config['nextseq']:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq'
        output:
            temp('Preprocessing/{type}.r1.trimmed.fq'),
            temp('Preprocessing/{type}.se1.trimmed.fq'),
            temp('Preprocessing/{type}.r2.trimmed.fq'),
            temp('Preprocessing/{type}.se2.trimmed.fq')
        threads:
            getThreads(10)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB,
        params:
            adapter = lambda wildcards: config["trimmomatic"]["adapter"][wildcards.type]
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing/trimming.{type}.log"
        message:
            "trimming: Trimmming {wildcards.type} reads."
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
            'Preprocessing/{type}.{read}.trimmed.fq'
        output:
            temp('Preprocessing/{type}.{read}.trimmed.fq'),
            temp('Preprocessing/{type}.{read}.trimmed.fq.gz')
        wildcard_constraints:
            read="se[1|2]"
        threads:
            getThreads(10)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing/nextseq_trimming_se.{type}.{read}.log"
        message:
            "trimming: Trimmming poly-G's off {wildcards.type} {wildcards.read} reads."
        shell:
            """
            cutadapt -j {threads} \
                     --nextseq-trim={config[trimmomatic][trailing]} \
                     -o {output[0]} \
                     --minimum-length {config[trimmomatic][minlen]} \
                     {input} >> {log} 2>&1
            pigz -p {threads} {input} >> {log} 2>&1
            """

    rule nextseq_trimming_pe:
        input:
            'Preprocessing/{type}.r1.trimmed.fq',
            'Preprocessing/{type}.r2.trimmed.fq'
        output:
            temp('Preprocessing/{type}.r1.trimmed.fq'),
            temp('Preprocessing/{type}.r2.trimmed.fq'),
            temp('Preprocessing/{type}.r1Cut.trimmed.fq'),
            temp('Preprocessing/{type}.r2Cut.trimmed.fq'),
            temp('Preprocessing/{type}.se1Cut.trimmed.fq'),
            temp('Preprocessing/{type}.se2Cut.trimmed.fq'),
            temp('Preprocessing/{type}.r1.trimmed.fq.gz'),
            temp('Preprocessing/{type}.r2.trimmed.fq.gz')
        threads:
            getThreads(10)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing/nextseq_trimming_pe.{type}.pe.log"
        message:
            "trimming: Trimmming poly-G's off {wildcards.type} pe reads."
        shell:
            """
            cutadapt -j {threads} \
                     --nextseq-trim={config[trimmomatic][trailing]} \
                     --minimum-length {config[trimmomatic][minlen]}:{config[trimmomatic][minlen]} \
                     -o {output[0]} \
                     -p {output[1]} \
                     --too-short-output {output[2]} \
                     --too-short-paired-output {output[3]} \
                     {input[0]} \
                     {input[1]} >> {log} 2>&1
            cutadapt -j {threads} \
                     --minimum-length {confitrimmo
                     --minimum-length {config[trimmomatic][minlen]} \
                     -o {output[5]} \
                     {output[3]} >> {log} 2>&1
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
            temp('Preprocessing/{type}.se.trimmed.fq'),
        threads:
            # Ensure at least 10gb of mem
            max(int(math.ceil(10 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing/cat_se_trimmed.{type}.log"
        message:
            "cat_se_trimmed: Concatenating trimmed single end {wildcards.type} reads"
        shell:
            """
            cat {input[0]} {input[1]} {input[2]} {input[3]} > {output} 2>> {log}
            pigz -p {threads} {input} >> {log} 2>&1
            """
else:
    rule trimming:
        input:
            'Preprocessing/{type}.r1.fq',
            'Preprocessing/{type}.r2.fq'
        output:
            temp('Preprocessing/{type}.r1.trimmed.fq'),
            temp('Preprocessing/{type}.se1.trimmed.fq'),
            temp('Preprocessing/{type}.r2.trimmed.fq'),
            temp('Preprocessing/{type}.se2.trimmed.fq')
        threads:
            getThreads(10)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB,
        params:
            adapter = lambda wildcards: config["trimmomatic"]["adapter"][wildcards.type],
            leading = config["trimmomatic"]["leading"],
            minlen = config["trimmomatic"]["minlen"],
            palindrome_clip_threshold = config["trimmomatic"]["palindrome_clip_threshold"],
            simple_clip_threshold = config["trimmomatic"]["simple_clip_threshold"],
            trailing = config["trimmomatic"]["trailing"],
            seed_mismatch = config["trimmomatic"]["seed_mismatch"],
            window_size = config["trimmomatic"]["window_size"],
            window_quality = config["trimmomatic"]["window_quality"],
            strictness = config["trimmomatic"]["strictness"],
            target_length = config["trimmomatic"]["target_length"]
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        log:
            "logs/preprocessing/trimming.{type}.log"
        message:
            "trimming: Trimmming {wildcards.type} reads."
        shell:
            """
            trimmomatic PE -threads {threads} \
                           {input[0]} \
                           {input[1]} \
                           {output} \
                           ILLUMINACLIP:{DBPATH}/adapters/{params.adapter}.fa:{params.seed_mismatch}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
                           LEADING:{params.leading} \
                           TRAILING:{params.trailing} \
                           SLIDINGWINDOW:{params.window_size}:{params.window_quality} \
                           MINLEN:{params.minlen} \
                           MAXINFO:{params.target_length}:{params.strictness} > {log} 2>&1
            """

    rule cat_se_trimmed:
        input:
            'Preprocessing/{type}.se1.trimmed.fq',
            'Preprocessing/{type}.se2.trimmed.fq'
        output:
            temp('Preprocessing/{type}.se.trimmed.fq')
        threads:
            # Ensure at least 10gb of mem
            max(int(math.ceil(10 / (int(MAIN_MEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]]))), 1)
        resources:
            runtime = SHORT_RUNTIME,
            mem_mb = MAIN_MEMCORE_MB
        log:
            "logs/preprocessing/cat_se_trimmed.{type}.log"
        message:
            "cat_se_trimmed: Concatenating trimmed single end {wildcards.type} reads"
        shell:
            """
            cat {input[0]} {input[1]} > {output} 2>> {log}
            """
