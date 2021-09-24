if config["preprocessing"]["nextseq"]:
    rule trimming:
        input:
            "Preprocessing/{omic}.r1.fq",
            "Preprocessing/{omic}.r2.fq"
        output:
            "Preprocessing/{omic}.r1.trimmoed.fq",
            "Preprocessing/{omic}.se1.trimmoed.fq",
            "Preprocessing/{omic}.r2.trimmoed.fq",
            "Preprocessing/{omic}.se2.trimmoed.fq"
        threads: 10
        resources:
            runtime="12:00:00",
            mem = get_mem_norm() * 10,
        params:
            params=get_trimmomatic_params
        conda:
            ENVDIR + "/IMP_preprocessing.yaml"
        log:
            "logs/preprocessing_trimming.{omic}.log"
        message:
            "trimming: Trimmming {wildcards.omic} reads."
        shell:
            "trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} {params.params} &> {log}"

    rule nextseq_trimming_se:
        input:
            "Preprocessing/{omic}.{read}.trimmoed.fq"
        output:
            "Preprocessing/{omic}.{read}.trimmed.fq",
        wildcard_constraints:
            read="se[1|2]"
        params:
            params="--nextseq-trim={} --minimum-length {}".format(
                config["preprocessing"]["trimmomatic"]["trailing"],
                config["preprocessing"]["trimmomatic"]["minlen"]
            ),
        threads: 10
        resources:
            runtime="12:00:00",
            mem = get_mem_norm() * 10,
        conda:
            ENVDIR + "/IMP_preprocessing.yaml"
        log:
            "logs/preprocessing_nextseq_trimming_se.{omic}.{read}.log"
        message:
            "trimming: Trimmming poly-Gs off {wildcards.omic} {wildcards.read} reads."
        shell:
            "cutadapt {params.params} -j {threads} -o {output} {input} &> {log}"

    rule nextseq_trimming_pe:
        input:
            "Preprocessing/{omic}.r1.trimmoed.fq",
            "Preprocessing/{omic}.r2.trimmoed.fq"
        output:
            "Preprocessing/{omic}.r1.trimmed.fq",
            "Preprocessing/{omic}.r2.trimmed.fq",
            "Preprocessing/{omic}.r1Cut.trimmed.fq",
            "Preprocessing/{omic}.r2Cut.trimmed.fq",
            "Preprocessing/{omic}.se1Cut.trimmed.fq",
            "Preprocessing/{omic}.se2Cut.trimmed.fq",
        threads: 10
        resources:
            runtime="12:00:00",
            mem = get_mem_norm() * 10,
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_nextseq_trimming_pe.{omic}.pe.log"
        message: "trimming: Trimmming poly-Gs off {wildcards.omic} pe reads."
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
            """
    rule cat_se_trimmed:
        input:
            "Preprocessing/{omic}.se1.trimmed.fq",
            "Preprocessing/{omic}.se2.trimmed.fq",
            "Preprocessing/{omic}.se1Cut.trimmed.fq",
            "Preprocessing/{omic}.se2Cut.trimmed.fq",
            "Preprocessing/{omic}.r1Cut.trimmed.fq",
            "Preprocessing/{omic}.r2Cut.trimmed.fq"
        output:
            "Preprocessing/{omic}.se.trimmed.fq",
        threads: 2
        resources:
            runtime="12:00:00",
            mem = get_mem_norm() * 2
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        message: "cat_se_trimmed: Concatenating trimmed single end {wildcards.omic} reads"
        shell:
            """
            cat {input[0]} {input[1]} {input[2]} {input[3]} > {output}
            pigz -p {threads} {input}
            """
else:
    rule trimming:
        input:
            "Preprocessing/{omic}.r1.fq",
            "Preprocessing/{omic}.r2.fq"
        output:
            "Preprocessing/{omic}.r1.trimmed.fq",
            "Preprocessing/{omic}.se1.trimmed.fq",
            "Preprocessing/{omic}.r2.trimmed.fq",
            "Preprocessing/{omic}.se2.trimmed.fq"
        threads: 10
        resources:
            runtime="12:00:00",
            mem = get_mem_norm() * 10,
        params:
            params=get_trimmomatic_params
        conda: ENVDIR + "/IMP_preprocessing.yaml"
        log: "logs/preprocessing_trimming.{omic}.log"
        message: "trimming: Trimmming {wildcards.omic} reads."
        shell:
            "trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} {params.params} &> {log}"

    rule cat_se_trimmed:
        input:
            "Preprocessing/{omic}.se1.trimmed.fq",
            "Preprocessing/{omic}.se2.trimmed.fq"
        output:
            "Preprocessing/{omic}.se.trimmed.fq"
        threads: 1
        resources:
            runtime="12:00:00",
            mem = get_mem_norm()
        message: "cat_se_trimmed: Concatenating trimmed single end {wildcards.omic} reads"
        shell:
            """
            cat {input[0]} {input[1]} > {output}
            """
