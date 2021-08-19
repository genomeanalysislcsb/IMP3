from os.path import getsize
KRAKSIZE = float(os.path.getsize( os.path.join(DBPATH, config['krakendb'], "hash.k2d") )) / 1000000000
KRAKTHREADS = round((KRAKSIZE / float(BIGMEMCORE[:-1])) + 0.5)

rule kraken_reads:
    input:
        r1="Preprocessing/{type}.r1.preprocessed.fq",
        r2="Preprocessing/{type}.r2.preprocessed.fq",
        se="Preprocessing/{type}.se.preprocessed.fq",
    output:
        labels_pe="Analysis/taxonomy/kraken/{type}.kraken.output.gz",
        labels_se="Analysis/taxonomy/kraken/{type}.se.kraken.output.gz",
        report_pe=report("Analysis/taxonomy/kraken/{type}.kraken.report", category="Taxonomy"),
        report_se=report("Analysis/taxonomy/kraken/{type}.se.kraken.report", category="Taxonomy"),
    log:
        "logs/analysis_kraken_reads.{type}.log"
    resources:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    params:
        # read_merger=os.path.join(SRCDIR, "read_merger.pl")
        db = os.path.join(DBPATH, config['krakendb']),
    threads:
        getThreads(KRAKTHREADS) 
    conda:
        os.path.join(ENVDIR,  "IMP_taxonomy.yaml")
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {config[krakendb]}."
    shell:
        """
        labels_pe={output.labels_pe}
        labels_se={output.labels_se}
        kraken2 -db {params.db} --threads {threads} --use-names --paired --output ${{labels_pe%.*}} --report {output.report_pe} {input.r1} {input.r2} &> {log}
        kraken2 -db {params.db} --threads {threads} --use-names --output ${{labels_se%.*}} --report {output.report_se} {input.se} &>> {log}
        gzip ${{labels_pe%.*}}
        gzip ${{labels_se%.*}}
        """
