# NOTE: based on https://git-r3lab.uni.lu/IMP/imp3/-/tree/issue38
# https://github.com/lmrodriguezr/nonpareil
#   does not support gzipped FASTQ files: https://github.com/lmrodriguezr/nonpareil/issues/35
rule preproc_nonpareil:
    input:
        "Preprocessing/mg.{reads}.preprocessed.fq"
    output:
        fq=temp("Preprocessing/mg.{reads}.preprocessed.nonpareil.fq"),
        out=expand("Stats/mg/mg.{{reads}}_nonpareil.{suffix}", suffix=["npa", "npc", "npl", "npo"])
    log:
        "logs/preproc_nonpareil.{reads}.log"
    benchmark:
        "benchmark/preproc_nonpareil.{reads}.tsv"
    wildcard_constraints:
        reads="|".join(get_read_types_raw()),
    params:
        minlen=config["preprocessing"]["trimmomatic"]["minlen"]
    threads: 1
    resources:
        runtime="08:00:00",
        mem=get_mem_norm,
    message:
        "Nonpareil: {input}"
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

rule preproc_nonpareil_curve:
    input:
        "Stats/mg/mg.{reads}_nonpareil.npo"
    output:
        report("Visualization/mg.{reads}_nonpareil.curve.png", category="Preprocessing"),
        report("Stats/mg/mg.{reads}_nonpareil.diversity.txt", category="Preprocessing")
    log:
        "logs/preproc_nonpareil_curve.{reads}.log"
    threads: 1
    resources:
        runtime="08:00:00",
        mem=get_mem_norm,
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    message:
        "Nonpareil: curve: {input}"
    script:
        os.path.join(SRCDIR, "nonpareil.R")
