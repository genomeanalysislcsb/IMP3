from os.path import getsize
KRAKSIZE = float(os.path.getsize( os.path.join(DBPATH, config['krakendb'], "hash.k2d") )) / 1000000000
KRAKTHREADS = round((KRAKSIZE / float(BIGMEMCORE[:-1])) + 0.5)

# https://ccb.jhu.edu/software/kraken2
# https://github.com/DerrickWood/kraken2/wiki/Manual
rule kraken_reads:
    input:
        "Preprocessing/{type}.se.preprocessed.fq"
    output:
        labels="Analysis/taxonomy/kraken/{type}.kraken.output.gz",
        report=report("Analysis/taxonomy/kraken/{type}.kraken.report", category="Taxonomy")
    resources:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    params:
        db = os.path.join(DBPATH, config['krakendb']),
    threads:
        getThreads(KRAKTHREADS) 
    conda:
        os.path.join(ENVDIR, "IMP_taxonomy.yaml")
    log:
        "logs/analysis_kraken_reads.{type}.log"
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {config[krakendb]}."
    shell:
        # --use-names: replace the tax. ID column with the scientific name and the tax. ID
        """
        labels={output.labels}
        kraken2 -db {params.db} --threads {threads} --use-names --output ${{labels%.*}} --report {output.report} {input} &> {log}
        gzip ${{labels%.*}}
        """