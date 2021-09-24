if config["preprocessing"]["nextseq"]:
    rule trimming:
        input:
            "Preprocessing/{omic}.se.fq"
        output:
            "Preprocessing/{omic}.se.trimmoed.fq"
        params:
            params=get_trimmomatic_params
        threads: 10
        resources:
            runtime = "12:00:00",
            mem = get_mem_norm() * 10,
        conda:
            ENVDIR + "/IMP_preprocessing.yaml"
        log:
            "logs/preprocessing_trimming.{omic}.log"
        message:
            "trimming: Trimmming {wildcards.omic} reads."
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.params} &> {log}"

    rule nextseq_trimming:
        input:
            "Preprocessing/{omic}.{read}.trimmoed.fq"
        output:
            "Preprocessing/{omic}.{read}.trimmed.fq",
        params:
            params="--nextseq-trim={} --minimum-length {}".format(
                config["preprocessing"]["trimmomatic"]["trailing"],
                config["preprocessing"]["trimmomatic"]["minlen"]
            ),
        threads: 10
        resources:
            runtime = "12:00:00",
            mem = get_mem_norm() * 10,
        conda:
            ENVDIR + "/IMP_preprocessing.yaml"
        log:
            "logs/preprocessing_nextseq_trimming.{omic}.{read}.log"
        message:
            "trimming: Trimmming poly-Gs off {wildcards.omic} {wildcards.read} reads."
        shell:
            "cutadapt {params.params} -j {threads} -o {output} {input} &> {log}"
else:
    rule trimming:
        input:
            "Preprocessing/{omic}.se.fq"
        output:
            "Preprocessing/{omic}.se.trimmed.fq"
        params:
            params=get_trimmomatic_params
        threads: 10
        resources:
            runtime = "12:00:00",
            mem = get_mem_norm() * 10,
        conda:
            ENVDIR + "/IMP_preprocessing.yaml"
        log:
            "logs/preprocessing_trimming.{omic}.log"
        message:
            "trimming: Trimmming {wildcards.omic} reads."
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.params} &> {log}"
