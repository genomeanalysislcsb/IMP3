# NOTE: based on code from XXX
# https://github.com/lmrodriguezr/nonpareil
#   does not support gzipped FASTQ files: https://github.com/lmrodriguezr/nonpareil/issues/35
rule nonpareil:
    input:
        'Preprocessing/mg.{read}.preprocessed.fq'
    output:
        fq=temp('Preprocessing/mg.{read}.preprocessed.nonpareil.fq'),
        out=expand(
            'Stats/mg/mg.{{read}}_nonpareil.{suffix}',
            suffix=["npa", "npc", "npl", "npo"]
        )
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    params:
        minlen = config['trimmomatic']['minlen']
    log:
        "logs/preprocessing_nonpareil.{read}.log"
    message:
        "nonpareil: processed {wildcards.read} fastq file."
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
        'Stats/mg/mg.{read}_nonpareil.npo'
    output:
        report('Visualization/mg.{read}_nonpareil.curve.png', category="Preprocessing"),
        report('Stats/mg/mg.{read}_nonpareil.diversity.txt', category="Preprocessing")
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    log:
        "logs/preprocessing_nonpareil_curve.{read}.log"
    message:
        "nonpareil_curve: processed {wildcards.read} data."
    script:
        os.path.join(SRCDIR, "nonpareil.R")
