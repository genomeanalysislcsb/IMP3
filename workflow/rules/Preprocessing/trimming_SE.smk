if config["preprocessing"]["nextseq"]:
    rule preproc_trim:
        input:
            "Preprocessing/{omic}.se.fq"
        output:
            "Preprocessing/{omic}.se.trimmoed.fq"
        log:
            "logs/preproc_trim.{omic}.log"
        benchmark:
            "benchmark/preproc_trim.{omic}.tsv"
        params:
            params=get_trimmomatic_params
        threads: 10
        resources:
            runtime="12:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message:
            "Trimming: {input}"
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.params} &> {log}"

    rule preproc_trim_nextseq:
        input:
            "Preprocessing/{omic}.{read}.trimmoed.fq"
        output:
            "Preprocessing/{omic}.{read}.trimmed.fq"
        log:
            "logs/preproc_trim_nextseq.{omic}.{read}.log"
        benchmark:
            "benchmark/preproc_trim_nextseq.{omic}.{read}.tsv"
        params:
            params=get_cutadapt_nextset_params,
        threads: 10
        resources:
            runtime="12:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message:
            "Nextseq trimming: {input}"
        shell:
            "cutadapt {params.params} -j {threads} -o {output} {input} &> {log}"
else:
    rule preproc_trim:
        input:
            "Preprocessing/{omic}.se.fq"
        output:
            "Preprocessing/{omic}.se.trimmed.fq"
        log:
            "logs/preproc_trim.{omic}.log"
        benchmark:
            "benchmark/preproc_trim.{omic}.tsv"
        params:
            params=get_trimmomatic_params
        threads: 10
        resources:
            runtime="12:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message:
            "Trimming: {input}"
        shell:
            "trimmomatic SE -threads {threads} {input} {output} {params.params} &> {log}"
