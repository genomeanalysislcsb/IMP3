from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = max(1, round((KRAKSIZE / (float(HIGH_MEM_BIGMEMCORE[:-1]) * to_gb_dict[config['cluster_specs']["mem_unit"]])) * 1.15))

rule kraken_reads:
    input:
        "Preprocessing/{type}.r1.preprocessed.fq",
        "Preprocessing/{type}.r2.preprocessed.fq",
    output:
        "Analysis/taxonomy/kraken/{type}.kraken.output.gz",
        report("Analysis/taxonomy/kraken/{type}.kraken.report",category="Taxonomy")
    resources:
        runtime = SHORT_RUNTIME,
        mem_mb = HIGH_MEM_BIGMEMCORE_MB
    params:
        tmp_out = lambda wildcards: "Analysis/taxonomy/kraken/" + wildcards.type + ".kraken.output",
        confidence = config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb']),
        report_minimizer = '--report-minimizer-data' if config['kraken']['report_minimizer'] else ''
    threads:
        getThreads(KRAKTHREADS, bigmem=True, high_mem_cluster=True)
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis_kraken_reads.{type}.log"
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {params.database}."
    shell:
        """
        mkdir -p Analysis/taxonomy/kraken
        kraken2 --db {params.database} \
                --confidence {params.confidence} \
                --threads {threads} \
                --paired \
                --use-names \
                {params.report_minimizer} \
                --output {params.tmp_out} \
                --report {output[1]} \
                {input} > {log} 2>&1
        gzip {params.tmp_out}
        """
