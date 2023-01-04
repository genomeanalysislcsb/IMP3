from os.path import getsize, join
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['kraken']['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = round((KRAKSIZE / (float(BIGMEMCORE[:-1]) * to_gb_dict[config["mem"]["mem_unit"]])) + 0.5)

rule kraken_contigs:
    input:
        "Assembly/%s.assembly.merged.fa" % ASS
    output:
        "Analysis/taxonomy/kraken/%s.kraken.contig.output" % ASS,
        "Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS
    resources:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    params:
        confidence = config['kraken']['confidence'],
        database = os.path.join(DBPATH, config['kraken']['krakendb'])
    threads:
        getThreads(KRAKTHREADS, bigmem=True)
    conda:
        ENVDIR + "/IMP_taxonomy.yaml"
    log:
        "logs/analysis_kraken_contigs.log"
    message:
        "kraken_contigs: Running kraken2 on contigs against {params.database}."
    shell:	
        """
        kraken2 --db {params.database} \
	            --confidence {params.confidence} \
                --threads {threads} \
                --use-names \
                --report-zero-counts \
                --output {output[0]} \
                --report {output[1]} \
                {input} > {log} 2>&1
        """

rule kraken_parse:
    input:
        "Analysis/taxonomy/kraken/%s.kraken.contig.output" % ASS,
        "Analysis/taxonomy/kraken/%s.kraken.contig.report" % ASS
    output:
        "Analysis/taxonomy/kraken/%s.kraken.parsed.tsv" % ASS
    threads:
        1
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    conda:
        ENVDIR + "/IMP_taxonomy.yaml"
    log:
        "logs/analysis_kraken_parse.log"
    message:
        "kraken_parse: Parsing kraken2 output for contigs."
    shell:
        """
        {SRCDIR}/krakenContigEntropy_dynamic.py -c {input[0]} \
                                                -e 0.1 \
                                                -o {output} \
                                                -r {input[1]} > {log} 2>&1
        """


