# if done: kraken to table: taxonomy per gene
# if done: mapping to table: depth of coverage per gene in each sample - normalize by total reads (MG and MT); reads per sample?
# if mapping done: summarize mapping: % samples with any coverage (MG and MT)s

# rule cluster2tab:
#    input:
#        "Catalogue/catalogue_cluster.tsv"
#    output:
#        "Catalogue/catalogue_genesPerRep.tsv"
#    threads: 1
#    resources:
#        runtime = "12:00:00",
#        mem = MEMCORE
#    params:
#        sample = samples
#    log: "logs/catalogue_cluster2tab.log"
#    message: 'catalogue_cluster2tab: checking number of genes represented by each representative.'
#    threads: 1
#    script:
#        SRCDIR + "/cluster2tab.py"

rule gff2tab:
    input:
        "Catalogue/catalogue.gff"
    output:
        "Catalogue/catalogue_anno.tsv"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        dbs = config['hmm_DBs'].split()
    log: "logs/catalogue_gff2tab.log"
    message: 'gff2tab: extracting gff for catalogue info.'
    threads: 1
    script:
        SRCDIR + "/gff2tab.py"

rule catalogue_tab:
    input:
        getCatIn
    output:
        "Catalogue/catalogue.tsv"
    resources:
        runtime = "12:00:00",
        mem = BIGMEMCORE
    params:
        types = TYPES
    threads: 1
    log: "logs/catalogue_tab.log"
    conda: ENVDIR + "/IMP_multi.yaml"
    message: "catalogue_tab: Merging results of catalogue."
    script:
        SRCDIR + "/multi_catalogue_merge.R"

