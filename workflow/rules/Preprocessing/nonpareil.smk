# NOTE: based on https://git-r3lab.uni.lu/IMP/imp3/-/tree/issue38
# https://github.com/lmrodriguezr/nonpareil
#   does not support gzipped FASTQ files: https://github.com/lmrodriguezr/nonpareil/issues/35
rule nonpareil:
    input:
        "Preprocessing/mg.{reads}.preprocessed.fq"
    output:
        fq=temp("Preprocessing/mg.{reads}.preprocessed.nonpareil.fq"),
        out=expand(
            "Stats/mg/mg.{{reads}}_nonpareil.{suffix}",
            suffix=["npa", "npc", "npl", "npo"]
        )
    resources:
        runtime="08:00:00",
        mem=get_mem_norm,
    params:
        minlen=config["preprocessing"]["trimmomatic"]["minlen"]
    log:
        "logs/preprocessing_nonpareil.{reads}.log"
    message:
        "nonpareil: processed {wildcards.reads} fastq file."
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    shell:
        """
        bname={output.out[0]}
        if [ {params.minlen} -ge 24 ]; then
            echo "reads >= 24"
            ln -s $(realpath {input}) {output.fq}
        else
            echo "rm reads < 24"
            seqkit seq {input} -j {threads} -m 24 | zcat > {output.fq}
        fi 2> {log}
        nonpareil -s {output.fq} -T kmer -f fastq -b ${{bname%.*}} &>> {log}
        """ 

rule nonpareil_curve:
    input:
        "Stats/mg/mg.{reads}_nonpareil.npo"
    output:
        report("Visualization/mg.{reads}_nonpareil.curve.png", category="Preprocessing"),
        report("Stats/mg/mg.{reads}_nonpareil.diversity.txt", category="Preprocessing")
    resources:
        runtime="08:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    log:
        "logs/preprocessing_nonpareil_curve.{reads}.log"
    message:
        "nonpareil_curve: processed {wildcards.reads} data."
    script:
        os.path.join(SRCDIR, "nonpareil.R")
