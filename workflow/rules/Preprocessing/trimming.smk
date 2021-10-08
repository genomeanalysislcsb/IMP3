if config["preprocessing"]["nextseq"]:
    rule preproc_trim:
        input:
            r1="Preprocessing/{omic}.r1.fq",
            r2="Preprocessing/{omic}.r2.fq",
        output:
            expand("Preprocessing/{{omic}}.{reads}.trimmoed.fq", reads=["r1", "se1", "r2", "se2"]) # NOTE: order matters
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
            "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output} {params.params} &> {log}"

    rule preproc_trim_nextseq_se:
        input:
            "Preprocessing/{omic}.{read}.trimmoed.fq"
        output:
            "Preprocessing/{omic}.{read}.trimmed.fq"
        log:
            "logs/preproc_trim_nextseq_se.{omic}.{read}.log"
        benchmark:
            "benchmark/preproc_trim_nextseq_se.{omic}.{read}.tsv"
        wildcard_constraints:
            read="se[1|2]"
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

    rule preproc_trim_nextseq_pe:
        input:
            r1="Preprocessing/{omic}.r1.trimmoed.fq",
            r2="Preprocessing/{omic}.r2.trimmoed.fq"
        output:
            r1="Preprocessing/{omic}.r1.trimmed.fq",
            r2="Preprocessing/{omic}.r2.trimmed.fq",
            r1c="Preprocessing/{omic}.r1Cut.trimmed.fq",
            r2c="Preprocessing/{omic}.r2Cut.trimmed.fq",
            se1c="Preprocessing/{omic}.se1Cut.trimmed.fq",
            se2c="Preprocessing/{omic}.se2Cut.trimmed.fq",
        log:
            "logs/preproc_trim_nextseq_pe.{omic}.pe.log"
        benchmark:
            "benchmark/preproc_trim_nextseq_pe.{omic}.pe.tsv"
        params:
            params=get_cutadapt_nextset_params,
            minlen=config["preprocessing"]["trimmomatic"]["minlen"],
        threads: 10
        resources:
            runtime="12:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message:
            "Nextseq trimming: {input}"
        shell:
            "cutadapt -j {threads} {params.params} -o {output.r1} -p {output.r1} "
            "--too-short-output {output.r1c} --too-short-paired-output {output.r2c} {input.r1} {input.r2} &> {log} && "
            "cutadapt -j {threads} --minimum-length {params.minlen} -o {output.se1c} {output.r1c} &>> {log} && "
            "cutadapt -j {threads} --minimum-length {params.minlen} -o {output.se2c} {output.r2c} &>> {log}"
    
    rule preproc_cat_trimmed_se:
        input:
            se=expand("Preprocessing/{{omic}}.{reads}.trimmed.fq", reads=["se1", "se2", "se1Cut", "se2Cut"]),
            pe=expand("Preprocessing/{{omic}}.{reads}.trimmed.fq", reads=["r1Cut", "r2Cut"]) # TODO: see gzip comment below
        output:
            "Preprocessing/{omic}.se.trimmed.fq",
        threads: 5
        resources:
            runtime="04:00:00",
            mem=get_mem_norm,
        conda:
            os.path.join(ENVDIR, "IMP_preprocessing.yaml")
        message:
            "Cat trimmed SE reads: {input}"
        shell:
            "cat {input.se} > {output} && pigz -p {threads} {input}" # TODO: gzip should be part of the archive rule and not done here
else:
    rule preproc_trim:
        input:
            r1="Preprocessing/{omic}.r1.fq",
            r2="Preprocessing/{omic}.r2.fq"
        output:
            expand("Preprocessing/{{omic}}.{reads}.trimmed.fq", reads=["r1", "se1", "r2", "se2"]) # NOTE: order matters
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
            "trimmomatic PE -threads {threads} {input.r1} {input.r2} {output} {params.params} &> {log}"

    rule preproc_cat_trimmed_se:
        input:
            se1="Preprocessing/{omic}.se1.trimmed.fq",
            se2="Preprocessing/{omic}.se2.trimmed.fq"
        output:
            "Preprocessing/{omic}.se.trimmed.fq"
        threads: 1
        resources:
            runtime="01:00:00",
            mem=get_mem_norm,
        message:
            "Cat trimmed SE reads: {input}"
        shell:
            "cat {input.se1} {input.se2} > {output}"
