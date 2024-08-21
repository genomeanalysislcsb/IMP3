if "preprocessing" in IMP_STEPS and int(config['fastp']['length_filtering']['length_required'])>=24:
    nonpareil_cmd = """
                    # Run nonpareil, touch empty output file if nonpareil fails
                    nonpareil -s {input} \
                              -T kmer \
                              -f fastq \
                              -b {params.out1}.{wildcards.read}{params.out2} >> {log} 2>&1 \
                      || touch {output} >> {log} 2>&1
                    """
else:
   nonpareil_cmd = """
                   export TMPDIR="{params.tmp_dir}"
                   TMP_FILE=$(mktemp "no_short_XXXXXX.fq")
                   seqkit seq {input} \
                              -j {threads} \
                              -m 24 >> $TMP_FILE 2>> {log}
                   
                   # Run nonpareil, touch empty output file if nonpareil fails
                   nonpareil -s $TMP_FILE \
                             -T kmer \
                             -f fastq \
                             -b {params.out1}.{wildcards.read}{params.out2} >> {log} 2>&1 \
                     || touch {output} >> {log} 2>&1
                   """

rule nonpareil:
    input:
        'Preprocessing/mg.{read}.preprocessed.fq'
    output:
        expand('Stats/mg/mg.{{read}}_nonpareil.{suffix}',suffix=["npa","npc","npl","npo"])
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        tmp_dir=TMPDIR,
        out1 = "Stats/mg/mg",
        out2 = "_nonpareil" 
    threads:
        1
    log:
        "logs/preprocessing/nonpareil.{read}.log"
    message:
        "nonpareil: processed {wildcards.read} fastq file."
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    shell:
        nonpareil_cmd

rule nonpareil_curve:
    input:
        'Stats/mg/mg.{read}_nonpareil.npo'
    output:
        report('Visualization/mg.{read}_nonpareil.curve.png',category="Preprocessing"),
        report('Stats/mg/mg.{read}_nonpareil.diversity.txt',category="Preprocessing")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "IMP_nonpareil.yaml")
    log:
        "logs/preprocessing/nonpareil_curve.{read}.log"
    message:
        "nonpareil_curve: processed {wildcards.read} data."
    script:
        SRCDIR + "/nonpareil.R"
