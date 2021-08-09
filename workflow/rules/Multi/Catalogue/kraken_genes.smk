from os.path import getsize
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['catalogue']['kraken']['db'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = round( (KRAKSIZE / float(BIGMEMCORE[:-1])) + 0.5)

rule kraken_cat:
    input:
        "Catalogue/catalogue.fasta" 
    output:
        "Catalogue/catalogue_kraken.output",
        "Catalogue/catalogue_kraken.report" 
    resources:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    threads: getThreads(KRAKTHREADS)
    conda: ENVDIR + "/IMP_taxonomy.yaml"
    log: "logs/analysis_kraken_cat"
    message: "kraken_cat: Running kraken2 on genes from catalog against {config[catalogue][kraken][db]}."
    shell:	
        """
        kraken2 -db {DBPATH}/{config[catalogue][kraken][db]} --threads {threads} --use-names --report-zero-counts \
         --output {output[0]} --report {output[1]} {input} >> {log} 2>&1
        """

rule kraken_parse_cat:
    input:
        "Catalogue/catalogue_kraken.output",
        "Catalogue/catalogue_kraken.report"
    output:
        "Catalogue/catalogue_kraken.parsed.tsv" 
    threads: 1
    resources:
        runtime = "2:00:00",
        mem = MEMCORE
    conda: ENVDIR + "/IMP_taxonomy.yaml"
    log: "logs/analysis_kraken_parse_cat.log"
    message: "kraken_parse_cat: Parsing kraken2 output for genes."
    shell:
        """
        {SRCDIR}/krakenContigEntropy_dynamic.py -c {input[0]} -e 0.1 -o {output} -r {input[1]} > {log} 2>&1
        """


