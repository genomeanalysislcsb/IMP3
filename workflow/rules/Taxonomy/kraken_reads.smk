from os.path import getsize
KRAKSIZE = float(os.path.getsize(DBPATH + "/" + config['krakendb'] + "/hash.k2d" )) / 1000000000
KRAKTHREADS = round( (KRAKSIZE / float(BIGMEMCORE[:-1])) + 0.5)

rule kraken_reads:
    input:
        r1="Preprocessing/{type}.r1.preprocessed.fq",
        r2="Preprocessing/{type}.r2.preprocessed.fq",
        se="Preprocessing/{type}.se.preprocessed.fq",
    output:
        reads=temp("Analysis/taxonomy/kraken/{type}.fq"),
        labels="Analysis/taxonomy/kraken/{type}.kraken.output.gz",
        report=report("Analysis/taxonomy/kraken/{type}.kraken.report",category="Taxonomy"),
    log:
        "logs/analysis_kraken_reads.{type}.log"
    resources:
        runtime = "8:00:00",
        mem = BIGMEMCORE
    params:
        read_merger=os.path.join(SRCDIR, "read_merger.pl")
    threads:
        getThreads(KRAKTHREADS) 
    conda:
        os.path.join(ENVDIR,  "IMP_taxonomy.yaml")
    message:
        "kraken_reads: Running kraken2 on {wildcards.type} reads against {config[krakendb]}."
    shell:
        # using PE and SE reads together as input: https://github.com/DerrickWood/kraken2/issues/214
        """
        {params.read_merger} --fq --output-format=legacy {input.r1} {input.r2} > {output.reads} &&
        cat {input.se} >> {output.reads} &&
        ofile={output.labels} &&
        kraken2 -db {DBPATH}/{config[krakendb]} --threads {threads} --use-names --output ${{ofile%.*}} --report {output.report} {output.reads} &> {log} &&
        gzip ${{ofile%.*}}
        """
