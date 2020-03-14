sum_stats_input = ["status/workflow.done"]
if not "preprocessing" in IMP_STEPS:
    sum_stats_input.append("status/preprocessing.analysed")
if CONTIGS:
    sum_stats_input.append("status/assembly.analysed")

localrules: ALL_stats 

gathered = [report("Stats/all_stats.tsv",category="Stats"),
            "Stats/all_stats.Rdata"]
if "binning" in IMP_STEPS:
    gathered.append(report("Stats/all_bin_stats.tsv",category="Binning"))

rule gather_stats:
    input:
        sum_stats_input
    output:
        gathered
    params:
        runtime = "12:00:00",
        mem = MEMCORE,
        steps = IMP_STEPS,
        ass = ASS,
        types = TYPES,
        dbs = config["hmm_DBs"]
    log: "logs/stats_gather_stats.log"
    message: "gather_stats: Summarizing stats."
    threads: 1
    conda: ENVDIR + "/IMP_binning.yaml"
    script:
        SRCDIR + "/summarize_all_stats.R"
         
rule ALL_stats:
    input:
        gathered
    output:
        touch("status/stats.done")

