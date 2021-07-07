if "preprocessing" in IMP_STEPS and int(config['trimmomatic']['minlen'])>=24:
    nonpareil_cmd = """
        nonpareil -s {input} -T kmer -f fastq -b {params.out1}.{wildcards.read}{params.out2} >> {log}
        """
else:
   nonpareil_cmd = """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "no_short_XXXXXX.fq")
        seqkit seq {input} -j {threads} -m 24 >> $TMP_FILE
        nonpareil -s $TMP_FILE -T kmer -f fastq -b {params.out1}.{wildcards.read}{params.out2} >> {log}
        """ 

rule nonpareil:
    input:
        'Preprocessing/mg.{read}.preprocessed.fq'
    output:
        expand('Stats/mg/mg.{{read}}_nonpareil.{suffix}',suffix=["npa","npc","npl","npo"])
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    params:
        out1 = "Stats/mg/mg",
        out2 = "_nonpareil" 
    threads: 1
    log: "logs/preprocessing_nonpareil.{read}.log"
    message: "nonpareil: processed {wildcards.read} fastq file."
    conda: ENVDIR + "/IMP_nonpareil.yaml"
    shell:
        nonpareil_cmd

rule nonpareil_curve:
    input:
        'Stats/mg/mg.{read}_nonpareil.npo'
    output:
        report('Visualization/mg.{read}_nonpareil.curve.png',category="Preprocessing"),
        report('Stats/mg/mg.{read}_nonpareil.diversity.txt',category="Preprocessing")
    resources:
        runtime = "8:00:00",
        mem = MEMCORE
    threads: 1
    conda: ENVDIR + "/IMP_nonpareil.yaml"
    log: "logs/preprocessing_nonpareil_curve.{read}.log"
    message: "nonpareil_curve: processed {wildcards.read} data."
    script:
        SRCDIR + "/nonpareil.R"
