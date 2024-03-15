
rule normalize_read_coverage_pe:
    input:
        mt_r1 = 'Preprocessing/mt.r1.preprocessed.fq',
        mt_r2 = 'Preprocessing/mt.r2.preprocessed.fq'
    output:
        mt_r1_norm = temp('Preprocessing/mt.r1.preprocessed.normalized.fq'),
        mt_r2_norm = temp('Preprocessing/mt.r2.preprocessed.normalized.fq')
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(128, bigmem=False)
    params:
        target_cov = config["read_coverage_normalization"]["target_cov"],
        min_cov = config["read_coverage_normalization"]["min_cov"]
    conda:
        ENVDIR + "/IMP_bbmap.yaml"
    log:
        "logs/assembly/preprocessing_norm_read_cov_pe.log"
    message:
        "normalize_read_coverage: Normalizing PE read coverage"
    shell:
        """
        bbnorm.sh in={input.mt_r1} \
                  in2={input.mt_r2} \
                  out={output.mt_r1_norm} \
                  out2={output.mt_r2_norm} \
                  target={params.target_cov} \
                  min={params.min_cov} \
                  threads={threads} > {log} 2>&1
        """

rule normalize_read_coverage_se:
    input:
        mt_se = 'Preprocessing/mt.se.preprocessed.fq'
    output:
        mt_se_norm = temp('Preprocessing/mt.se.preprocessed.normalized.fq')
    resources:
        runtime = MEDIUM_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        getThreads(128, bigmem=False)
    params:
        target_cov = config["read_coverage_normalization"]["target_cov"],
        min_cov = config["read_coverage_normalization"]["min_cov"]
    conda:
        ENVDIR + "/IMP_bbmap.yaml"
    log:
        "logs/assembly/preprocessing_norm_read_cov_se.log"
    message:
        "normalize_read_coverage: Normalizing SE read coverage"
    shell:
        """
        bbnorm.sh in={input.mt_se} \
                  out={output.mt_se_norm} \
                  target={params.target_cov} \
                  min={params.min_cov} \
                  threads={threads} > {log} 2>&1
        """