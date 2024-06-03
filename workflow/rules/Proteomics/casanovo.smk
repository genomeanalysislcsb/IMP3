rule casanovo:
    input:
        config = "<some_path>/1229_2_sipros_test_all_samp_db_01/config/1229_2_sipros.cfg",
        ms_data = directory("<some_path>/1229_2_sipros_test_01/data"),
        db = rules.sipros_prepare_db.output.db
    output:
        out_dir = directory("<some_path>/1229_2_sipros_test_all_samp_db_01/out")
    resources:
        runtime = LONG_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    threads:
        1
    conda:
        os.path.join(ENVDIR, "casanovo.yaml")
    log:
        "logs/analysis/run_sipros.log"
    message:
        "Running sipros."
    shell:
        """
        casanovo sequence -o results.mztab spectra.mgf > {log} 2>&1
        """