rule collate_stats:
    input:
        expand("{raw_directory}/Stats/all_stats.tsv", raw_directory=samples.path[samples.stats == 1])
    output:
        report("Collection/multi_stats.tsv",category="Collection"),
        "Collection/multi_stats.RDS"
    threads: 1
    resources:
        runtime = "12:00:00",
        mem = MEMCORE
    params:
        samples = samples.loc[samples.stats == 1,"sample"]
    log: "logs/collate_stats.log"
    message: 'collate_stats: stats from ' + " ".join(samples.loc[samples.stats == 1,"sample"])+ '.'
    threads: 1
    conda: ENVDIR + "/IMP_multi.yaml"
    script:
        SRCDIR + "/multi_collect_stats.R"


localrules:  ctrl_stats_collation

rule ctrl_stats_collation:
    input:
        "Collection/multi_stats.tsv"
    output:
        touch("status/stats.collated")


