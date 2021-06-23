def getFeats(wildcards):
    all=expand("{raw_directory}/Analysis/annotation/{type}.{db}.counts.tsv",raw_directory=samples.loc[[a and b for a,b in zip(samples[wildcards.db] == 1, samples[wildcards.type] == 1)],"path"].tolist(),type=wildcards.type,db=wildcards.db)
    return all

rule collate_featureCounts:
    input:
        getFeats
    output:
        "Collection/multi_{type}.{db}.tsv",
        "Collection/multi_{type}.{db}.RDS"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        samples = lambda wildcards: samples.loc[[a and b for a,b in zip(samples[wildcards.db] == 1, samples[wildcards.type] == 1)],"sample"].tolist()
    log: "logs/collate_featureCounts_{type}.{db}.log"
    message: 'collate_featureCounts: collecting {wildcards.db} {wildcards.type} data.'
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/multi_collect_featureCounts.R"

localrules:  ctrl_featureCounts_collation

rule ctrl_featureCounts_collation:
    input:
        expand("Collection/multi_{type}.{db}.tsv",db=config['hmm_DBs'].split(),type=TYPES)
    output:
        touch("status/features.collated")
