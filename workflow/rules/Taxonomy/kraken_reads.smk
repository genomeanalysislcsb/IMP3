from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = max(1, round((KRAKSIZE / (float(HIGH_MEM_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]])) * 1.15))

rule kraken_reads:
    input:
        "Preprocessing/{type}.r1.preprocessed.fq",
        "Preprocessing/{type}.r2.preprocessed.fq",
    output:
        "Analysis/taxonomy/kraken/{type}.pe.kraken.output.gz",
        report("Analysis/taxonomy/kraken/{type}.pe.kraken.report", category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = HIGH_MEM_BIGMEMCORE_MB
    params:
        tmp_out = lambda wildcards: f"Analysis/taxonomy/kraken/{wildcards.type}.pe.kraken.output",
        confidence = config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb']),
        report_minimizer = '--report-minimizer-data' if config['kraken']['report_minimizer'] else '',
        min_hit_groups = config['kraken']['min_hit_groups']
    threads:
        getThreads(KRAKTHREADS, bigmem=True, high_mem_cluster=True)
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis/taxonomy/kraken_reads.{type}.pe.log"
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {params.database}."
    shell:
        """
        mkdir -p Analysis/taxonomy/kraken >> {log} 2>&1
        kraken2 --db {params.database} \
                --confidence {params.confidence} \
                --threads {threads} \
                --paired \
                --use-names \
                --minimum-hit-groups {params.min_hit_groups} \
                {params.report_minimizer} \
                --output {params.tmp_out} \
                --report {output[1]} \
                {input} >> {log} 2>&1
        gzip {params.tmp_out} >> {log} 2>&1
        """

rule combine_kraken_reports:
    input:
        pe_report = "Analysis/taxonomy/kraken/{type}.pe.kraken.report",
        se_report = "Analysis/taxonomy/kraken/{type}.se.kraken.report"
    output:
        report("Analysis/taxonomy/kraken/{type}.kraken.report",category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = MAIN_MEMCORE_MB
    params:
        script_path = os.path.join(ROOTDIR, "submodules", "KrakenTools", "combine_kreports.py"),
        combined_report_options = "--only-combined --no-headers"
    log:
        "logs/analysis/taxonomy/combine_kraken_reports.{type}.log"
    message:
        "Combining Kraken reports for {wildcards.type} into a single report."
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    threads:
        getThreads(12, bigmem=False, high_mem_cluster=False)
    shell:
        """
        python {params.script_path} \
            -r {input.se_report} {input.pe_report} \
            -o {output} \
            {params.combined_report_options} > {log} 2>&1
        """

