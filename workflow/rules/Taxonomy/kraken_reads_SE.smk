from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = round((KRAKSIZE / (float(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]])) + 0.5)


rule kraken_reads:
    input:
        "Preprocessing/{type}.se.preprocessed.fq"
    output:
        "Analysis/taxonomy/kraken/{type}.kraken.output.gz",
        report("Analysis/taxonomy/kraken/{type}.kraken.report",category="Taxonomy")
    resources:
        runtime = "08:00:00",
        mem_mb = BIGMEMCORE_MB
    params:
        tmp_out= lambda wildcards: "Analysis/taxonomy/kraken/" + wildcards.type + ".kraken.output",
        confidence=config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb'])
    threads:
        getThreads(KRAKTHREADS, bigmem=True)
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
	            --use-names \
	            --output {params.tmp_out} \
	            --report {output[1]} {input} > {log} 2>&1
	    gzip {params.tmp_out}
        """

